import sys
sys.path.append("../preprocess")
from qbiolip import load_db, aa_3to1
import argparse
import yaml

arg_parser = argparse.ArgumentParser(description="""
Build inputs files (.yaml) for Boltz-2 predictor.
Use this script after using "build_fasta.py" (MSA should be pre-computed first).
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
    default="../../databases/boltz2/inputs"
)
arg_parser.add_argument(
    "--msa-path",
    help="Path to pre-computed MSA using build_fasta.py file. Default: ../../databases/boltz2/msa",
    default="/work/lepikhov/ml_docking_review/databases/boltz2/msa"
)
a = arg_parser.parse_args()
class CustomDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(CustomDumper, self).increase_indent(flow, False)

df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
df = df.loc[df.target]

with open("./all_fastas.fasta", "r") as fasta_f:
    all_fasta = [l.replace("\n", "") for l in fasta_f]
msa_names = [l.replace(">", "") for l in all_fasta if l.startswith(">")]
sequences = [l for l in all_fasta if not l.startswith(">")]
seq_msa = {seq:msa for msa, seq in zip(msa_names, sequences)} # mapping of sequence - corresponding MSA file
for i in range(df.shape[0]):
    row = df.iloc[i]
    for is_trimed in [True, False]:
        rec_path = a.rec_trimed_path if is_trimed else a.rec_path
        assembly_filename = row.ligand_detail if is_trimed else row.assembly_id
        with open(f"{rec_path}/{assembly_filename}.pdb") as rec_f:
            res_data = [(
                aa_3to1[r[17:20].strip()], #residue
                r[21] # chain
            ) for r in rec_f if r[:6].strip() == "ATOM" and r[12:16].strip() == "CA"]
        chain_seq_dict = {chain:"" for _,chain in res_data} # initialize all chains as empty sequences
        for res,chain in res_data: # fill the chain's sequence
            chain_seq_dict[chain] = chain_seq_dict[chain] + res
        name = row.ligand_detail + ("", "_trimed")[is_trimed]
        data = {
            "sequences": [
                {   
                    "protein": {
                        "id": chain,
                        "sequence": seq,
                        "msa": f"{a.msa_path}/{seq_msa[seq]}.a3m",
                    }
                }
            for chain, seq in chain_seq_dict.items()]
        }
        data["sequences"].append({
            "ligand": {
                "id": "B",
                "smiles": row.lig_canon_smiles
            }
        })
        yaml_str = yaml.dump(data, Dumper=CustomDumper, sort_keys=False, indent=2, default_flow_style=False)
        out_path = f"{a.output_path}/{name}.yaml"
        with open(out_path, "w") as out_f:
            out_f.write(yaml_str)
