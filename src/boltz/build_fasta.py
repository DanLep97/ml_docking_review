import sys
sys.path.append("../preprocess")
from qbiolip import load_db, aa_3to1
import argparse

arg_parser = argparse.ArgumentParser(description="""
Build sequence inputs for Boltz-2 predictor for MSA pre-computing.
Use this file to compute MSA.
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


all_fastas = ""
processed_sequences = []
for i in range(df.shape[0]):
    row = df.iloc[i]
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
    for sequence, chain in u_chain_seq_dict.items():
        if sequence in processed_sequences:
            continue
        processed_sequences.append(sequence)
        all_fastas += f">{row.assembly_id}-{chain}\n"
        all_fastas += sequence + "\n"

with open("./all_fastas.txt", "w") as fasta_f:
    fasta_f.write(all_fastas)