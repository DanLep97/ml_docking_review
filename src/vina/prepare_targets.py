import sys
sys.path.append("../preprocess")
from qbiolip import load_db
import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from vina import Vina
import torch
import os
import glob
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

arg_parser = argparse.ArgumentParser(
    description="""
    Prepares targets to generate docking results using vina.
    To run this script, first install AutDockTools on your computer.
    """
)
arg_parser.add_argument(
    "--qbiolip",
    help="Path to qbiolip database. Default: ../../databases/QBioLiP/qbiolip.csv",
    default="../../databases/QBioLiP/qbiolip.csv"
)
arg_parser.add_argument(
    "--ec_uniprot_ids",
    help="Path to ec_uniprot_ids.txt file. Default: ../../databases/ec_uniprot_ids.txt",
    default="../../databases/ec_uniprot_ids.txt"
)
arg_parser.add_argument(
    "--rec-path",
    help="Path to receptors. Default: ../../databases/QBioLiP/nr_rec",
    default="../../databases/QBioLiP/nr_rec"
)
arg_parser.add_argument(
    "--rec-trimed-path",
    help="Path to trimed receptors. Default: ../../databases/QBioLiP/nr_rec_trimed",
    default="../../databases/QBioLiP/nr_rec_trimed"
)
arg_parser.add_argument(
    "--lig-path",
    help="Path to nr_lig. Default: ../../databases/QBioLiP/nr_lig",
    default="../../databases/QBioLiP/nr_lig"
)
arg_parser.add_argument(
    "--output",
    help="Path to docked structures using Vina. Default: ../../databases/Vina/outputs",
    default="../../databases/Vina/outputs"
)
arg_parser.add_argument(
    "--adfr",
    help="""Path to ADFR bin folder.
    Download from:  
    https://ccsb.scripps.edu/adfr/downloads/
    This package is used to prepare ligands and receptors for Vina (adding H + converting to pdbqt).
    """,
    default="/home/daniil/Desktop/Ph.D/ADFR/bin"
)
a = arg_parser.parse_args()

df = load_db(a.qbiolip, a.ec_uniprot_ids)
df = df.loc[df.target]
df.reset_index(inplace=True, drop=True)

banned_cases = ["3ooj_2", "4bbj_1", "8avy_1", "4gz5_5"]
# subprocess.run(f"rm -r {a.output}/*", shell=True)
print(f"Cleared old outputs from {a.output}")
subprocess.run(f"rm {a.rec_path}/*.pdbqt", shell=True)
subprocess.run(f"rm {a.rec_trimed_path}/*.pdbqt", shell=True)
subprocess.run(f"rm {a.lig_path}/*.pdbqt", shell=True)
print("Removed pdbqt files.")
failed_pdbqt = 0
failed_docking = 0
failed_pbd_docking = 0
failed_embeddings = 0
curr_pdqt = True # use this variable to go to next case
curr_path = os.getcwd()

for i in range(df.shape[0]):
# for i in range(3):
    case = df.loc[i] 
    name = f"{a.output}/{case.ligand_detail}"
    if os.path.exists(name) and os.path.exists(f"{name}_trimed") or case.assembly_id in banned_cases:
        continue

    # add H to receptors using pdbfixer:
    # print("Adding H to receptors (trimed and whole) and converting to pdqt...")
    receptor_paths = [
        f"{a.rec_path}/{case.assembly_id}.pdb",
        f"{a.rec_trimed_path}/{case.ligand_detail}.pdb"
    ]
    for rec_p in receptor_paths:
        receptor_H = rec_p.replace('.pdb', '_reduced.pdb')
        receptor_pdbqt = receptor_H.replace('_reduced.pdb', ".pdbqt")
        # Add H using reduce then convert to pdbqt file:
        if not os.path.exists(receptor_pdbqt):
            print(f"Converting {rec_p.split('/')[-1]} file to pdbqt...")
            # Add H:
            subprocess.run(
                f"{a.adfr}/reduce -Quiet -DB {a.adfr}/reduce_wwPDB_het_dict.txt {rec_p} > {receptor_H}",
                shell=True
            )
            # Convert to pdbqt:
            subprocess.run(
                f"{a.adfr}/prepare_receptor -r {receptor_H} -o {receptor_pdbqt}",
                shell=True
            )
            # remove intermediate file:
            subprocess.run(f"rm {receptor_H}", shell=True)
        print("Build receptor pdbqt:", receptor_pdbqt.split("/")[-1])
 
    #add H to ligands using adfr package:
    print("Adding H to ligands and converting to pdbqt..")
    ligand_path = f"{case.ligand_detail}_randconf.pdb"
    ligand_pdbqt = ligand_path.replace(".pdb", ".pdbqt")
    os.chdir(a.lig_path)
    subprocess.run(f"{a.adfr}/prepare_ligand -l {ligand_path} -o {ligand_pdbqt} -A hydrogens", shell=True)
    os.chdir(curr_path)
            
    # dock structure: 
    # from vina paper, it is recommended that the size of the box is size of ligand + 15A
    # or higher than 22.5 total.
    # Here, since we're attempting blind docking, vina grid size is fixed and is:
    # max length on each axis + 20A.
    # The center is the COM of the protein or of the ligand (TBD or PBD).
    with open(f"{a.lig_path}/{ligand_path}") as lig_f:
        lig_lines = [l for l in lig_f if l[:6].strip() == "HETATM"]
    lig_xyz = torch.tensor([
        [
            float(r[30:38].strip()), #x
            float(r[38:46].strip()), #y
            float(r[46:54].strip())  #z
        ]
    for r in lig_lines])
    lig_xyz_expd = lig_xyz.expand(lig_xyz.shape[0], lig_xyz.shape[0], -1)
    # maximum distances between atoms, per axis:
    axis_dist = (lig_xyz_expd - lig_xyz_expd.transpose(0,1)).abs().max(dim=0).values.max(dim=0).values
    # maximum size among all 3 axis:
    max_axis_size = axis_dist.max()
    for rec_p in receptor_paths:
        receptor_pdbqt = rec_p.replace(".pdb", ".pdbqt")
        name = case.ligand_detail if "_trimed" not in rec_p else case.ligand_detail + "_trimed"
        out_dir = f"{a.output}/{name}"

        if "_trimed" not in name: # TBD setup, use receptor coordinates to build grid
            with open(rec_p) as rec_f:
                rec_lines = [l for l in rec_f if l[:6].strip() == "ATOM" and l[76:78].strip() != "H"]
            rec_xyz = torch.tensor([
                [
                    float(r[30:38].strip()), #x
                    float(r[38:46].strip()), #y
                    float(r[46:54].strip())  #z
                ]
            for r in rec_lines])
            centroid_xyz = rec_xyz.mean(dim=0).tolist() # the center of the grid is the center of the protein
            # compute the size of the the box around the whole protein:
            rec_xyz_expd = rec_xyz.expand(rec_xyz.shape[0], rec_xyz.shape[0], -1)
            axis_dist = (rec_xyz_expd - rec_xyz_expd.transpose(0,1)).abs().max(dim=0).values.max(dim=0).values
            box_size = (axis_dist + max_axis_size).float().tolist()
            print(f"BUILDING TBD with {name}")
            
        else: # PBD setup, use ligand coordinates to build grid
            centroid_xyz = lig_xyz.mean(dim=0).tolist()
            # the box size is this maximum distance repeated for each axis:
            box_size = (max_axis_size.tile(3)+max_axis_size/2).tolist()
            print(f"BUILDING PBD with {name}")

        print("Box size:", box_size)
        print("Centroid:", centroid_xyz)
        if len(glob.glob(f"{out_dir}/*")) > 1:# or not curr_pdqt:
            continue
        try:
            os.mkdir(out_dir)
        except:
            print(f"{out_dir} already created, skipping folder creation.")
        
        try:
            v = Vina(sf_name="vina")
            v.set_receptor(receptor_pdbqt)
            v.set_ligand_from_file(f"{a.lig_path}/{ligand_pdbqt}")
            v.compute_vina_maps(
                center=centroid_xyz,
                box_size=box_size
            )
            # score the pose:
            energy = v.score()
            print(f"Score before energy minimization: {energy}")
            
            # dock the ligand:
            v.dock(exhaustiveness=32, n_poses=20)
            if "_trimed" in rec_p:
                out_f = f"{out_dir}/{case.ligand_detail}_trimed.pdbqt"
            else:
                out_f = f"{out_dir}/{case.ligand_detail}.pdbqt"
            v.write_poses(out_f, n_poses=20, overwrite=True)
            return_code = subprocess.run(f"obabel {out_f} -O {out_f.replace('.pdbqt', '.pdb')} -m", shell=True)
        except Exception as e:
            if "_trimed" not in name:
                failed_docking += 1
            else:
                failed_pbd_docking += 1
            print(f"failed docking. Error: {e}")
print("Job done.")
print("Number of failed docking: ", failed_docking)
print("Number of failed PBD docking: ", failed_pbd_docking)