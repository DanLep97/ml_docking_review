import sys
sys.path.append("../preprocess")
from qbiolip import load_db
import argparse
import os
import subprocess
import pandas

arg_parser = argparse.ArgumentParser(description="""
    Prepare targets for DynamicBind.
    https://github.com/luwei0917/DynamicBind/tree/main
""")
arg_parser.add_argument(
    "--qbiolip",
    help="Path to qbiolip database. Default: ../../databases/QBioLiP/qbiolip.csv",
    default="../../databases/QBioLiP/qbiolip.csv"
)
arg_parser.add_argument(
    "--ec-uniprot-ids",
    help="Path to ec_uniprot_ids.txt file. Default: ../../databases/ec_uniprot_ids.txt",
    default="../../databases/ec_uniprot_ids.txt"
)
arg_parser.add_argument(
    "--nr-rec",
    help="Path to nr_rec directory. Default: ../../databases/QBioLiP/nr_rec",
    default="/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_rec"
)
arg_parser.add_argument(
    "--nr-rec-trimed",
    help="Path to nr_rec_trimed directory. Default: ../../databases/QBioLiP/nr_rec_trimed",
    default="/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_rec_trimed"
)
arg_parser.add_argument(
    "--nr-lig",
    help="Path to nr_lig directory. Default: ../../databases/QBioLiP/nr_lig",
    default="/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_lig"
)
a = arg_parser.parse_args()

df = load_db(a.qbiolip, a.ec_uniprot_ids)
df = df.loc[df.target]
df.reset_index(drop=True, inplace=True)

protein_paths = []
ligand_paths = []
names = []
for i in range(df.shape[0]):
    case = df.iloc[i]
    protein_path = [
        f"{a.nr_rec}/{case.assembly_id}.pdb",
        f"{a.nr_rec_trimed}/{case.ligand_detail}.pdb",
    ]
    name = [
        f"{case.ligand_detail}",
        f"{case.ligand_detail}_trimed"
    ]
    ligand_path = f"{a.nr_lig}/{case.ligand_detail}_randconf.sdf"
    if not os.path.exists(ligand_path):
        pdb2sdf_command_code = subprocess.run(
            f"obabel {ligand_path.replace('.sdf', '.pdb')} -O {ligand_path}",
        shell=True, check=True)
    ligand_path = [ligand_path]*2
    protein_paths.extend(protein_path)
    ligand_paths.extend(ligand_path)
    names.extend(name)

out_df = pandas.DataFrame({
    "name": names,
    "ligand": ligand_paths,
    "protein_path": protein_paths
})

cases_per_files = 10
for j,i in enumerate(range(0,out_df.shape[0], cases_per_files)):
    start = j*cases_per_files
    end = start+cases_per_files
    curr_df = out_df.iloc[start:end].to_csv(f"./csv/ligands_{j}.csv", index=False)
print(f"Saved protein_ligand.csv with {out_df.shape[0]} entries.")