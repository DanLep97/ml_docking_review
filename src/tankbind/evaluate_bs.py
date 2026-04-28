import sys
sys.path.append("../preprocess")
from qbiolip import load_db
import argparse
import pandas
import torch

arg_parser = argparse.ArgumentParser(description="""
Use this script to evaluate the accuracy of predicted
binding site in comparison between TankBind and p2rank.
A correctly predicted center of a binding site (hit) is 
considered if the center of the predicted BS is in a 
cutoff distance (default: 5A).

To use this script, provide the output of tankbind's csv
(the info.csv file).
""")

arg_parser.add_argument(
    "--qbiolip",
    help="Path to qbiolip db.",
    default="../../databases/QBioLiP/qbiolip.csv"
)
arg_parser.add_argument(
    "--uniprot-ids",
    help="Path to ec_uniprot_ids.txt file. Default: ../../databases/ec_uniprot_ids.txt",
    default="../../databases/ec_uniprot_ids.txt"
)
arg_parser.add_argument(
    "--info",
    help="The csv containing the prediction of tankbind",
    default="../../databases/TankBind/outputs/info_with_BA.csv"
)
arg_parser.add_argument(
    "--nr-lig",
    help="Path to nr_rec folder. Default: ../../databases/QBioLiP/nr_rec",
    default="../../databases/QBioLiP/nr_lig"
)
a = arg_parser.parse_args()

db = load_db(a.qbiolip, a.uniprot_ids)
df = db[db.target]
print(df.ligand_detail)
info = pandas.read_csv(a.info)

top1_dists = []
com1_dists = []
for i in range(df.shape[0]):
    try:
        # retrieve top1 and com1
        case = df.iloc[i]
        info_row = info.loc[info.protein_name == case.ligand_detail]
        com1 = info_row.loc[info_row.pocket_name == "pocket_1"].pocket_com.values[0].split(",")
        com1 = torch.tensor([float(c) for c in com1])
        top1 = info_row.sort_values(by=["affinity"], ascending=False).iloc[0].pocket_com.split(",")
        top1 = torch.tensor([float(c) for c in top1])

        # retrieve center of ligand:
        with open(f"{a.nr_lig}/{case.ligand_detail}.pdb") as lig_f:
            xyz_lig = torch.tensor([
                [
                    float(r[30:38].strip()), #x
                    float(r[38:46].strip()), #y
                    float(r[46:54].strip())  #z
                ] for r in lig_f if r[0:6].strip() == "HETATM"
            ])
            com = xyz_lig.mean(dim=0)
        # com = torch.tensor([
        #     case.centroid_x,
        #     case.centroid_y,
        #     case.centroid_z
        # ])

        top1_dists.append((top1 - com).pow(2).sum().sqrt())
        com1_dists.append((com1 - com).pow(2).sum().sqrt())
    except Exception as e:
        continue

top1_dists = torch.stack(top1_dists)
com1_dists = torch.stack(com1_dists)
print("top1 dists shape:", top1_dists.shape)
print("com1 dists shape:", com1_dists.shape)
print("5 top1:", top1_dists[0:5])

print(f"Hits for top1 TankBind prediction: {(top1_dists < 2).nonzero().shape[0]/top1_dists.shape[0]*100:.2f}%")
print(f"Hits for com1 p2rank prediction: {(com1_dists < 2).nonzero().shape[0]/com1_dists.shape[0]*100:.2f}%")
# print("Avg distance for com1:", com1_dists.mean())
# print(top1_dists)
