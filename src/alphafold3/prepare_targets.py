import sys
sys.path.append("../preprocess")
from qbiolip import load_db, aa_3to1
import json
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import glob

arg_parser = argparse.ArgumentParser(description="""
Builds the .json file containing inputs for alphafold3.
Use this file to compute MSA for single chains that later will be combined.
""")
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
    "--output-path",
    help="Path where the json files will be saved.",
    default="../../databases/alphafold3/inputs"
)
a = arg_parser.parse_args()

df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
df = df.loc[df.target]

skip_cases = [
    "8avy_1_Y01_J".upper(), "8avy_1_Y01_J_trimed".upper(),
    "4ayx_1_ACP_D1".upper(), "4ayx_1_ACP_D1_trimed".upper(),
    "5och_4_Y01_BA".upper(), "5och_4_Y01_BA_trimed".upper()
]

all_sequences = []
all_cases = 0
all_names = []
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
    # get sequence for untrimed receptor:
    with open(f"{a.rec_path}/{row.assembly_id}.pdb") as rec_f:
        res_data = [(
            aa_3to1[r[17:20].strip()], #residue
            r[21] # chain
        ) for r in rec_f if r[:6].strip() == "ATOM" and r[12:16].strip() == "CA"]
    chain_seq_dict = {chain:"" for _,chain in res_data} # initialize all chains as empty sequences
    for res,chain in res_data: # fill the chain's sequence
        chain_seq_dict[chain] = chain_seq_dict[chain] + res
    u_chain_seq_dict = {v:k for k,v in chain_seq_dict.items()}
    name = row.ligand_detail
    all_names.append(name)
    modelSeeds = [0]
    sequences = [{
        "protein": {
            "id": chain,
            "sequence": sequence,
        } 
    } for sequence, chain in u_chain_seq_dict.items() if sequence not in all_sequences]
    if len(sequences) == 0 or any(len(s["protein"]["sequence"]) > 900 for s in sequences):
        continue
    all_sequences.extend([s["protein"]["sequence"] for s in sequences])

    entry = {
        "name": name,
        "modelSeeds": modelSeeds,
        "sequences": sequences,
        "dialect": "alphafold3",
        "version": 1
    }
    if name.upper() not in skip_cases:
        json.dump(entry, open(f"{a.output_path}/{name}.json","w"))
        all_cases += 1
print("all cases: ", all_cases)