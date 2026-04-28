import sys
sys.path.append("../preprocess")
from qbiolip import load_db, aa_3to1
import json
import argparse
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from repo import pdb_to_ccd_mmcif

arg_parser = argparse.ArgumentParser(description="""
Builds the .json file containing inputs for alphafold3.
""")
arg_parser.add_argument(
    "--lig-path",
    help="Path to the ligand folder. Default: ../../databases/QBioLiP/nr_lig",
    default="../../databases/QBioLiP/nr_lig",
)
arg_parser.add_argument(
    "--rec-path",
    help="Path to the receptor folder. Default: ../../databases/QBioLiP/nr_rec",
    default="../../databases/QBioLiP/nr_rec",
)
arg_parser.add_argument(
    "--rec-trimed-path",
    help="Path to the receptor folder. Default: ../../databases/QBioLiP/nr_rec_trimed",
    default="../../databases/QBioLiP/nr_rec_trimed",
)
arg_parser.add_argument(
    "--input-path",
    help="Path where the json of unique sequences are saved.",
    default="../../databases/alphafold3/unique_chain_outputs"
)
arg_parser.add_argument(
    "--af-output",
    help="Path to alphafold3 output file.",
    default="../../databases/alphafold3/multi_chain_outputs"
)
arg_parser.add_argument(
    "--cif-path",
    help="Path to custom CCD MMCIF files (ligand format with precise xyz coordinates).",
    default="../../databases/alphafold3/cif"
)
a = arg_parser.parse_args()

df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
df = df.loc[df.target]

unique_sequences_files = glob.glob(f"{a.input_path}/*/*.json")

u_sequences = []
msas = []

for i,f_path in enumerate(unique_sequences_files):
    # if i ==0:
    with open(f_path) as f:
        entry = json.load(f)["sequences"]
    u_sequences.extend(s["protein"]["sequence"] for s in entry)
    msas.extend(s["protein"]["unpairedMsa"] for s in entry)
print("n u_sequences:", len(u_sequences))
print("n Msas:", len(msas))
sequence_msa = {s:m for s,m in zip(u_sequences,msas)}
modelSeeds = list(range(1))
all_cases=0
print("n u_sequences in dict:", len(sequence_msa))

for i in range(df.shape[0]):
    row = df.iloc[i]
    smiles = row.lig_canon_smiles
    # check if ligand can be loaded:
    options = AllChem.ETKDGv3()
    options.clearConfs = False
    conf_id = -1
    try:
        mol = Chem.MolFromSmiles(row.lig_canon_smiles)
        mol = Chem.AddHs(mol)
        canonical_order = AllChem.CanonicalRankAtoms(mol)
        for atom, can_idx in zip(mol.GetAtoms(), canonical_order):
            atom.SetProp("name", atom.GetSymbol().upper() + str(can_idx + 1))
        conf_id = AllChem.EmbedMolecule(mol, options)
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
        mol = Chem.AddHs(mol)
        if mol is None:
            continue
        if conf_id != -1:
            conformer = mol.GetConformer(conf_id)
    except:
        continue
    for trimed in [True, False]:
        rec_path = f"{a.rec_path}/{row.assembly_id}.pdb" if not trimed else f"{a.rec_trimed_path}/{row.ligand_detail}.pdb"
        with open(rec_path) as rec_f:
            res_data = [(
                aa_3to1[r[17:20].strip()], #residue
                r[21] # chain
            ) for r in rec_f if r[:6].strip() == "ATOM" and r[12:16].strip() == "CA"]
        chain_seq_dict = {chain:"" for _,chain in res_data} # initialize all chains as empty sequences
        for res,chain in res_data: # fill the chain's sequence
            chain_seq_dict[chain] = chain_seq_dict[chain] + res
        name = row.ligand_detail + ("","_trimed")[trimed]
        sequences = [{
            "protein": {
                "id": chain,
                "sequence": sequence,
                "unpairedMsa": sequence_msa[sequence],
                "pairedMsa": "",
                "templates": [],
            }, 
        } for chain, sequence in chain_seq_dict.items() if sequence in u_sequences]
        seq_len =sum([len(s["protein"]["sequence"]) for s in sequences])
        if seq_len > 1400:
            continue
        ligand_id = [i for i in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if i not in list(chain_seq_dict.keys())][0]
        # convert pdb to ccd mmcif:
        lig_path = f"{a.lig_path}/{row.ligand_detail}_randconf.pdb"
        mccif_file = f"{a.cif_path}/{name}.cif"
        pdb_to_ccd_mmcif(lig_path, mccif_file, "MY_LIG")
        with open(mccif_file) as mf:
            mccif_output = "\n".join([l for l in mf])

        sequences.append({
            "ligand": {
                "id": ligand_id,
                "ccdCodes": ["MY_LIG"]
            }
        })

        entry = dict(
            name=name,
            modelSeeds = modelSeeds,
            sequences = sequences,
            dialect="alphafold3",
            version=1,
            userCCD=mccif_output
        )

        out_path = f"{a.af_output}/{name}.json"
        if not os.path.exists(f"../../databases/alphafold3/outputs/{name.lower()}"):
            json.dump(entry, open(out_path,"w"), indent=4)
            all_cases += 1
print("number of cases:", all_cases)
