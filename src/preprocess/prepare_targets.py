if __name__ != "__main__":
    from .qbiolip import load_db
else:
    from qbiolip import load_db
import argparse
import pandas
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import torch
import subprocess
aa_3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    "SEC": "U"
}

aa_1to3 = {value:key for key,value in aa_3to1.items()}

def savePdb(original_pdb_file_name, new_coordinates, new_pdb_name='result.pdb'):
    write_str = ""
    for i, line in enumerate(open(original_pdb_file_name).readlines()):
        if line[:6].strip() == "HETATM":
            newline =    line[:30] + ('%.3f' % new_coordinates[i][0]).rjust(8) + line[38:]
            newline = newline[:38] + ('%.3f' % new_coordinates[i][1]).rjust(8) + newline[46:]
            newline = newline[:46] + ('%.3f' % new_coordinates[i][2]).rjust(8) + newline[54:]
            write_str = write_str + newline 
    with open(new_pdb_name, "w") as out:
        out.write(write_str)

def rotation_matrix(alpha, beta, gamma):
    # Create rotation matrices for each axis
    R_x = torch.tensor([
        [1, 0, 0],
        [0, torch.cos(alpha), -torch.sin(alpha)],
        [0, torch.sin(alpha), torch.cos(alpha)]
    ])

    R_y = torch.tensor([
        [torch.cos(beta), 0, torch.sin(beta)],
        [0, 1, 0],
        [-torch.sin(beta), 0, torch.cos(beta)]
    ])

    R_z = torch.tensor([
        [torch.cos(gamma), -torch.sin(gamma), 0],
        [torch.sin(gamma), torch.cos(gamma), 0],
        [0, 0, 1]
    ])

    # Combine the rotations: R = R_z * R_y * R_x
    R = torch.mm(R_z, torch.mm(R_y, R_x))

    return R

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description="""
        Given the qbiolip database, annotates enzyme-ligand complexes as targets.
        Targets are used to evaluate docking softwares and tools.
        An enzyme-ligand complex is considered as a target if the receptor file exists
        and if the ligand can be loaded with rdkit.
        """,
    )
    arg_parser.add_argument(
        "--ntargets",
        help="""
        Targets are enzyme-ligand complexes used as inputs for docking tools and softwares.
        This arguments specifies how many complexes are selected per class. Default: 30
        """,
        default=60,
        type=int
    )
    arg_parser.add_argument(
        "--cutoff",
        help="""
        Distance cutoff (in Angstrom) from binding residues after which all residues are removed.
        Default: 10.
        """,
        default=10,
        type=int
    )
    arg_parser.add_argument(
        "--lig-min-atoms",
        help="""
        Number of minimum atoms for the ligand when selecting targets.
        Default: 5
        """,
        default=5,
        type=int
    )
    arg_parser.add_argument(
        "--rec-min-res", 
        help="""
        Number of mininum target's enzyme residues. Default: 100
        """,
        default=100,
        type=int
    )
    arg_parser.add_argument(
        "--rec-max-res", 
        help="""
        Number of maximum target's enzyme residues. Default: 1600
        """,
        default=1600,
        type=int
    )
    arg_parser.add_argument(
        "--ec-uniprot-ids",
        help="""
        Path to biolip csv. Default: ../../databases/ec_uniprot_ids.txt 
        """,
        default="../../databases/ec_uniprot_ids.txt"
    )
    arg_parser.add_argument(
        "--csv-path",
        help="""
        Path to the annotation file. Default: ../../databases/QBioLiP/PL_nr.csv
        """,
        default="../../databases/QBioLiP/qbiolip.csv" 
    )

    a = arg_parser.parse_args()
    alpha = torch.randint(180, (1,)).deg2rad()
    beta = torch.randint(180, (1,)).deg2rad()
    gamma = torch.randint(180, (1,)).deg2rad()
    R = rotation_matrix(alpha, beta, gamma)

    df = load_db(a.csv_path, a.ec_uniprot_ids, no_metal_ligands=True)
    df["target"] = False
    groups = df.group.dropna().unique()
    lig_folder = "../../databases/QBioLiP/nr_lig"
    rec_folder = "../../databases/QBioLiP/nr_rec"

    # first, groupby assembly id and ligand id
    df_unique_ligands = df.groupby(["assembly_id","ligand_id"]).filter(lambda x: len(x) == 1)
    # then, only by assembly id so each assembly ID is modelled only once.
    # This line is not useful for now, keeping track of the ids in the pdbs should be enough
    # df_unique_assembly_id = df_unique_ligands.groupby(["assembly_id"]).filter(lambda x: len(x) == 1)

    # Filter out elements from the equibind pdbs:
    prohibited_or_already_picked_pdbs = [l.replace("\n", "") for l in open("../../databases/EquiBind/equibind_pdbs.txt")]
    
    subprocess.run(f"rm {lig_folder}/*randconf.pdb", shell=True)
    subprocess.run(f"rm {lig_folder}/*randconf.sdf", shell=True)
    target_indices = []
    print("groups: ",groups)
    for group in groups:
        group_indices = pandas.Series(df.loc[
            (df.index.isin(df_unique_ligands.index)) &
            (df.group == group) &
            (df.lig_atoms.str.len() > a.lig_min_atoms) & 
            (df.rec_len_resi < a.rec_max_res) & 
            (df.rec_len_resi > a.rec_min_res) &
            (
                df.stoichiometry.str.contains("Hetero") | 
                df.stoichiometry.str.contains("Monomer")
            ) &
            (~df.lig_atoms.isna())
        ].index).sample(frac=1).tolist()
        n = 0 #number of built targets per group

        while n < a.ntargets:
            # verify the ligand:
            try:
                ind = group_indices.pop(0)
            except:
                print(f"Not enough indices left. n: {n}, group: {group}")
                exit()
            case = df.loc[ind]
            if case.pdb_id in prohibited_or_already_picked_pdbs:
                continue
            lig_detail = case.ligand_detail
            ligand_src = f"{lig_folder}/{lig_detail}.pdb"
            receptor_src = f"{rec_folder}/{case.assembly_id}.pdb"
            receptor_trimed_dest = f"{rec_folder}_trimed/{lig_detail}.pdb"

            try:
                if os.path.isfile(ligand_src) and os.path.isfile(receptor_src):
                    with open(ligand_src) as lig_f:
                        lig_lines = [l for l in lig_f]
                        elements = set([r[76:78].strip() for r in lig_lines if r[:6].strip() == "HETATM"])
                        if "H" in elements:
                            continue # work only with ligand without H
                        # check that it is not a polymer: atom names are not found more than ounce:
                        target_atom_names = [
                            r[12:16].strip() for r in lig_lines
                        ]
                        if len(target_atom_names) != len(set(target_atom_names)):
                            continue # repeating names of the same atom = potential polymer.
                            # it's also a problem to have same atom name for AF3
                    lig = Chem.MolFromPDBFile(ligand_src, removeHs=True)
                    AllChem.EmbedMolecule(lig)
                    AllChem.UFFOptimizeMolecule(lig)
                    lig_out = f"{lig_folder}/{lig_detail}_test.pdb"
                    with open(lig_out, "w") as f:
                        f.write(Chem.MolToPDBBlock(lig))
                    # check if the ligand pdb can be generated, if not, remove the case:
                    lig_check = Chem.MolFromPDBFile(lig_out, sanitize=True, removeHs=False)
                    if lig_check is None:
                        os.remove(lig_out)
                        continue
                    else:
                        lig_check = None
                        n=n+1
                        target_indices.append(ind)
                        prohibited_or_already_picked_pdbs.append(case.pdb_id)
                else:
                    continue
            except:
                continue

            # trim the receptor from all aa further from a cutoff distance:
            with open(receptor_src) as rec_f:
                # get the xyz of all alpha carbon atoms:
                rec_atoms = [l.replace("\n", "") for l in rec_f if l[:6].strip() == "ATOM"]
            ca = [l for l in rec_atoms if l[12:16].strip() == "CA"]
            resi = [int(l[22:26].strip()) for l in ca]
            rec_xyz = [
                [
                    float(r[30:38].strip()), #x
                    float(r[38:46].strip()), #y
                    float(r[46:54].strip())  #z
                ] for r in ca
            ]
            rec_xyz = torch.tensor(rec_xyz)
            rec_chains = pandas.Series([l[21] for l in ca]) # filtering by index is easier with a series
            # get the xyz coordinates of the atoms:
            with open(ligand_src) as lig_f:
                # get the xyz of all ligand atoms:
                lig_atoms = [l.replace("\n", "") for l in lig_f if l[:6].strip() == "HETATM"]
            lig_xyz = torch.tensor([
                [
                    float(l[30:38].strip()), #x
                    float(l[38:46].strip()), #y
                    float(l[46:54].strip()), #z
                ] for l in lig_atoms
            ])

            lig_m = lig_xyz.expand(rec_xyz.shape[0], lig_xyz.shape[0], -1)
            rec_m = rec_xyz.expand(lig_xyz.shape[0], rec_xyz.shape[0], -1)
            # pairwise distance matrix:
            dist = (lig_m - rec_m.transpose(0,1)).norm(dim=2)
            in_cutoff = (dist < a.cutoff).nonzero()
            # we are interested only in alpha carbons which are in the cutoff distance from the ligand.
            # given the shape of the pairwise distance matrix, the index in the first dimension is the 
            # index of the atom in the list of atoms:
            in_cutoff_atoms = in_cutoff[:,0].unique() # indices of the atoms
            in_cutoff_chains = rec_chains[in_cutoff_atoms.tolist()].unique() # list of chains in the cutoff distance
            new_atoms = [record for record in rec_atoms if record[21] in in_cutoff_chains]
            with open(receptor_trimed_dest, "w") as trimed_f:
                trimed_f.write("\n".join(new_atoms))
            # define the ligand center of mass:
            ligand_com = lig_xyz.mean(dim=0)
            # define rotation matrix:
            alpha = torch.randint(180, (1,)).deg2rad()
            beta = torch.randint(180, (1,)).deg2rad()
            gamma = torch.randint(180, (1,)).deg2rad()
            r = rotation_matrix(alpha, beta, gamma)
            # center ligand xyz coordinates, apply random rotation and translate back to its original center of mass:
            rand_lig_xyz = (lig_xyz - ligand_com)@r + ligand_com
            # modify the ligand xyz coordinates:
            savePdb(lig_out, rand_lig_xyz, f"{lig_folder}/{lig_detail}_randconf.pdb")

            #recalculate the new size of grid:
            size = rand_lig_xyz.max(dim=0).values
            size[size < 22.5] = 22.5

            df.loc[ind, ["size_x", "size_y", "size_z"]] = [
                size[0].item(), size[1].item(), size[2].item()
            ]
            # print("Path to receptor:", receptor_src)
            # print("Path to trimed receptor:", receptor_trimed_dest)
            # print("Path to the ligand: ", ligand_src)
            # print("shape of all ca xyz:", rec_xyz.shape)
            # print("shape of all ligand atoms:", lig_xyz.shape)
            # print("shape of expanded lig_xyz: ", lig_m.shape)
            # print("shape of expanded rec_xyz: ", rec_m.shape)
            # print("shape of dist: ", dist.shape)
            # print("len of all alpha carbons inside the cutoff:", in_cutoff_atoms.shape[0])
            # print("alpha carbons inside the cutoff:", in_cutoff_atoms)

    # select targets:
    df.loc[target_indices, ["target"]] = True
    df.to_csv(a.csv_path, index=False)
    print(df.loc[df.target])
    print(f"Database {a.csv_path} updated with {df.loc[df.target].shape[0]} targets ({a.ntargets} times {groups.shape[0]}).")
    print(f"Last target index:", target_indices[-1])
    print("Unique assembly IDs:", df[df.target].assembly_id.unique().shape[0])