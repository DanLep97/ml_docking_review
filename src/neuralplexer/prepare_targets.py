import pandas
import argparse
import sys
sys.path.append("../preprocess")
from qbiolip import load_db
import glob

arg_parser = argparse.ArgumentParser(
    description="""
    Prepare targets for NeuralPLexer.
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
    help="Path to NeuralPLexer predictions. Default: ../../databases/NeuralPLexer/inferences",
    default="../../databases/NeuralPLexer/inferences"
)
a = arg_parser.parse_args()
df = load_db(a.qbiolip, a.ec_uniprot_ids)
df = df.loc[df.target]
df.reset_index(drop=True, inplace=True)

rec_paths = []
lig_paths = []
output_paths = []
modeled_paths = glob.glob(f"{a.output_path}/*")
already_modeled = 0
for i in range(df.shape[0]):
    case = df.loc[i]
    lig_path = f"{a.lig_input}/{case.ligand_detail}_randconf.pdb"
    output_path = f"{a.output_path}/{case.ligand_detail}"
    output_path_trimed = f"{a.output_path}/{case.ligand_detail}_trimed" 

    # add the true blind receptor and ligand:
    if output_path not in modeled_paths:
        rec_paths.append(f"{a.rec_input}/{case.assembly_id}.pdb")
        lig_paths.append(lig_path)
        output_paths.append(output_path)
    else:
        already_modeled+=1

    # add the pseudo blind receptor and ligand:
    if output_path_trimed not in modeled_paths:
        rec_paths.append(f"{a.trimed_rec_input}/{case.ligand_detail}.pdb")
        lig_paths.append(lig_path)
        output_paths.append(output_path_trimed)
    else:
        already_modeled+=1
print("Already modeled:", already_modeled)

with open("rec_paths.txt", "w") as rec_f:
    rec_f.write("\n".join(rec_paths))
with open("lig_paths.txt", "w") as lig_f:
    lig_f.write("\n".join(lig_paths))
with open("output_paths.txt", "w") as output_f:
    output_f.write("\n".join(output_paths))