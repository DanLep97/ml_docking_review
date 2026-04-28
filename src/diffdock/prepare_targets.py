import sys
sys.path.append("../preprocess")
from qbiolip import load_db
import pandas
import subprocess
import argparse

arg_parser = argparse.ArgumentParser(
    description="""
    Prepare targets for DiffDock.
"""
)
arg_parser.add_argument(
    "--lig-src",
    help="""
    Path to ligand source. Default: ../../databases/QBioLiP/nr_lig
    """,
    default="/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_lig"
)
arg_parser.add_argument(
    "--rec-src",
    help="""
    Path to receptor source. Default: ../../databases/QBioLiP/nr_rec
    """,
    default="/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_rec"
)
arg_parser.add_argument(
    "--nrows",
    help="Number of rows in the target csv file. Default: 10",
    default=10,
    type=int
)
arg_parser.add_argument(
    "--diffdock-csv",
    help="""
    Path to diffdock csv file used by diffdock for inference.
    Default: /home/daniil/Desktop/Ph.D/DiffDock/data/targets.csv
    """,
    default="/home/daniil/Desktop/Ph.D/DiffDock/data/targets"
)
a = arg_parser.parse_args()

df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
df = df.loc[df.target]
df.reset_index(drop=True,inplace=True)
indices = df.index.tolist()

# columns for diffdock df:
complex_names = []
protein_paths = []
ligand_descriptions = []
protein_sequence = None

for ind in indices:
    case = df.loc[ind]
    ligand_src = f"{a.lig_src}/{case.ligand_detail}_randconf.pdb"
    ligand_out = f"{a.lig_src}/{case.ligand_detail}_randconf.sdf"

    pdb2sdf_command_code = subprocess.run(
        f"obabel {ligand_src} -O {ligand_out}",
    shell=True, check=True)

    # add whole protein lines to the csv:
    complex_names.append(case.ligand_detail)
    protein_paths.append(f"{a.rec_src}/{case.assembly_id}.pdb")
    ligand_descriptions.append(ligand_out)

    # add trimed protein lines to the csv:
    complex_names.append(f"{case.ligand_detail}_trimed")
    protein_paths.append(f"{a.rec_src}_trimed/{case.ligand_detail}.pdb")
    ligand_descriptions.append(ligand_out)

diffdock_df = pandas.DataFrame({
    "complex_name": complex_names,
    "protein_path": protein_paths,
    "ligand_description": ligand_descriptions,
    "protein_sequence": protein_sequence
})

nsplits = int(diffdock_df.shape[0]/a.nrows)
for split in range(nsplits):
    start = a.nrows*split
    end = start+a.nrows
    if end > diffdock_df.shape[0]:
        end = None
    diffdock_df[start:end].to_csv(f"{a.diffdock_csv}_{split}.csv", index=False)