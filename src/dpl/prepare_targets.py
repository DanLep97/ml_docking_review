import pandas
import argparse
import sys
sys.path.append("../preprocess")
from qbiolip import load_db
import os
import subprocess
import glob

arg_parser = argparse.ArgumentParser(
    description="""
    Prepare targets for DPL.
"""
)
arg_parser.add_argument(
    "--qbiolip",
    help="Path to the qbiolip database. Default: ../../databases/QBioLiP/qbiolip.csv",
    default="../../databases/QBioLiP/qbiolip.csv",
)
arg_parser.add_argument(
    "--ec-uniprot-ids",
    help="Path to the ec_uniprot_ids.txt. Default: ../../databases/ec_uniprot_ids.txt",
    default="../../databases/ec_uniprot_ids.txt"
)
arg_parser.add_argument(
    "--rec-input",
    help="Path to nr_rec. Default: ../../databases/QBioLiP/nr_rec",
    default="../../databases/QBioLiP/nr_rec"
)
arg_parser.add_argument(
    "--trimed-rec-input",
    help="Path to nr_rec_trimed. Default: ../../databases/QBioLiP/nr_rec_trimed",
    default="../../databases/QBioLiP/nr_rec_trimed"
)
arg_parser.add_argument(
    "--lig-input",
    help="Path to nr_lig. Default: ../../databases/QBioLiP/nr_lig",
    default="../../databases/QBioLiP/nr_lig"
)
arg_parser.add_argument(
    "--output-path",
    help="Path to NeuralPLexer predictions. Default: ../../databases/DPL/inferences",
    default="../../databases/DPL/outputs/"
)
a = arg_parser.parse_args()
df = load_db(a.qbiolip, a.ec_uniprot_ids)
df = df.loc[df.target]
df.reset_index(drop=True, inplace=True)
existing_paths = glob.glob(f"{a.output_path}/*")
existing_cases = [p.split("/")[-1] for p in existing_paths]
modeled_cases = 0

rec_paths = []
lig_paths = []
output_paths = []
for i in range(df.shape[0]):
    case = df.loc[i]
    if case.ligand_detail in existing_cases:
        modeled_cases+=1
        continue
    lig_input = f"{a.lig_input}/{case.ligand_detail}_randconf.sdf"
    if not os.path.exists(lig_input):
        pdb2sdf_command_code = subprocess.run(
            f"obabel {lig_input.replace('.sdf', '.pdb')} -O {lig_input}",
        shell=True, check=True)
    # add the true blind receptor and ligand:
    rec_paths.append(f"{a.rec_input}/{case.assembly_id}.pdb")
    lig_paths.append(lig_input)
    output_paths.append(f"{a.output_path}/{case.ligand_detail}")

    # add the pseudo blind receptor and ligand:
    rec_paths.append(f"{a.trimed_rec_input}/{case.ligand_detail}.pdb")
    lig_paths.append(lig_input)
    output_paths.append(f"{a.output_path}/{case.ligand_detail}_trimed")

with open("rec_paths.txt", "w") as rec_f:
    rec_f.write("\n".join(rec_paths))
with open("lig_paths.txt", "w") as lig_f:
    lig_f.write("\n".join(lig_paths))
with open("output_paths.txt", "w") as output_f:
    output_f.write("\n".join(output_paths))
print("Total number of modelled cases:", modeled_cases)