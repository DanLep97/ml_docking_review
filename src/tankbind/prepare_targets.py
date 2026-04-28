import sys
from Bio.PDB import PDBParser
sys.path.append("../")
from preprocess.qbiolip import load_db
from preprocess.prepare_targets import savePdb
from rdkit import Chem
import argparse
import pandas as pd
import os
import subprocess
import torch
import logging
from torch_geometric.loader import DataLoader
from tqdm import tqdm
torch.multiprocessing.set_sharing_strategy('file_system')

arg_parser = argparse.ArgumentParser(
    description="""
    Preprocess data and generate results for TankBind.
"""
)
arg_parser.add_argument(
    "--tankbind",
    help="Path to tankbind package. Default: /home/daniil/Desktop/Ph.D/TankBind/tankbind",
    default="/home/daniil/Desktop/Ph.D/TankBind/tankbind"
)
arg_parser.add_argument(
    "--tankbind-model",
    help="Path to self_dock.pt model. Default: /home/daniil/Desktop/Ph.D/TankBind/saved_models/self_dock.pt",
    default="/home/daniil/Desktop/Ph.D/TankBind/saved_models/self_dock.pt"
)
arg_parser.add_argument(
    "--qbiolip",
    help="Path to qbiolip database. Default: ../../databases/QBioLiP/qbiolip.csv",
    default="../../databases/QBioLiP/qbiolip.csv"
)
arg_parser.add_argument(
    "--uniprot-ids",
    help="Path to ec_uniprot_ids.txt file. Default: ../../databases/ec_uniprot_ids.txt",
    default="../../databases/ec_uniprot_ids.txt"
)
arg_parser.add_argument(
    "--nr-rec",
    help="Path to nr_rec folder. Default: ../../databases/QBioLiP/nr_rec",
    default="../../databases/QBioLiP/nr_rec"
)
arg_parser.add_argument(
    "--nr-rec-trimed",
    help="Path to nr_rec_trimed folder. Default: ../../databases/QBioLiP/nr_rec_trimed",
    default="../../databases/QBioLiP/nr_rec_trimed"
)
arg_parser.add_argument(
    "--nr-lig",
    help="Path to nr_lig folder. Default: ../../databases/QBioLiP/nr_lig",
    default="../../databases/QBioLiP/nr_lig"
)
arg_parser.add_argument(
    "--output-path",
    help="Output path. Default: ../../databases/TankBind/outputs",
    default="../../databases/TankBind/outputs"
)
arg_parser.add_argument(
    "--ligand-output-path",
    help="Path where the computed ligand's xyz coordinates are saved. Default: ../../databases/TankBind/outputs/lig_pdb",
    default="../../databases/TankBind/outputs/lig_pdb"
)
arg_parser.add_argument(
    "--dataset-path",
    help="Path for the TankBind database to saved at. Default: ../../databases/TankBind/dataset",
    default="../../databases/TankBind/dataset"
)
arg_parser.add_argument(
    "--p2rank",
    help="Path to prank binary. Default: /home/daniil/Downloads/p2rank_2.3/prank",
    default="/home/daniil/Downloads/p2rank_2.5/prank"
)
arg_parser.add_argument(
    "--p2rank-out",
    help="Path to the output of the p2rank binary. Default: ../../databases/TankBind/p2rank",
    default="../../databases/TankBind/p2rank"
)
a = arg_parser.parse_args()
sys.path.insert(0, a.tankbind)
# tankbind imports:
from feature_utils import get_protein_feature
from feature_utils import extract_torchdrug_feature_from_mol
from data import TankBind_prediction
from model import get_model
from generation_utils import get_LAS_distance_constraint_mask, get_info_pred_distance

df = load_db(a.qbiolip,a.uniprot_ids)
df = df.loc[df.target]
df.reset_index(drop=True, inplace=True)
# get binding residue info
rec_list = [
    [
        f"{a.nr_rec}/{df.iloc[n].assembly_id}.pdb",
        f"{a.nr_rec_trimed}/{df.iloc[n].ligand_detail}.pdb"
    ] for n in range(df.shape[0])
]
rec_list = set([p for n in rec_list for p in n])
with open("./prank.txt", "w") as p_f:
    p_f.write("\n".join(rec_list))
br_cmd = f"{a.p2rank} predict ./prank.txt -o {a.p2rank_out} -threads 32"
os.system(br_cmd)
print("Binding site predictions ready.")
print("Building protein and ligand features..")
protein_dict = {}
compound_dict = {}
all_smiles = []
pdbs = []
failed_cases = 0
for i in range(df.shape[0]):
    # generate receptor feature:
    case = df.iloc[i]
    protein_files = [f"{a.nr_rec}/{case.assembly_id}.pdb", f"{a.nr_rec_trimed}/{case.ligand_detail}.pdb"]
    for protein_file in protein_files:
        ligand_file = f"{a.nr_lig}/{case.ligand_detail}_randconf.sdf"
        parser = PDBParser(QUIET=True)
        protein_key = case.ligand_detail
        if "_trimed" in protein_file:
            protein_key = protein_key + "_trimed"
        s = parser.get_structure(protein_key, protein_file)
        res_list = list(s.get_residues())
        try:
            protein_dict[protein_key] = get_protein_feature(res_list)
        except Exception as e:
            failed_cases+=1
            continue

        # generate ligand feature:
        if not os.path.exists(ligand_file):
            pdb2sdf_command_code = subprocess.run(
                f"obabel {ligand_file.replace('.sdf', '.pdb')} -O {ligand_file}",
            shell=True, check=True)

        mol = Chem.MolFromMolFile(ligand_file)
        try:
            compound_dict[protein_key] = extract_torchdrug_feature_from_mol(mol, has_LAS_mask=True)
        except Exception as e:
            del protein_dict[protein_key]
            failed_cases+=1
            continue
        all_smiles.append(case.lig_canon_smiles)
        pdbs.append(protein_file.split("/")[-1].replace(".pdb",""))
print("Failed cases:", failed_cases)

# build database:
print("Building dataset...")
info = []
rec_lig_zip = zip(protein_dict.keys(), compound_dict.keys())
for i, (pdb, compound_name) in enumerate(rec_lig_zip):
    # use protein center as the block center.
    com = ",".join([str(a.round(3)) for a in protein_dict[pdb][0].mean(axis=0).numpy()])
    # info.append([pdb, compound_name, "protein_center", com])

    p2rankFile = f"{a.p2rank_out}/{pdbs[i]}.pdb_predictions.csv"
    pocket = pd.read_csv(p2rankFile)
    pocket.columns = pocket.columns.str.strip()
    pocket_coms = pocket[['center_x', 'center_y', 'center_z']].values
    for ith_pocket, com in enumerate(pocket_coms):
        # if ith_pocket == 0:
        com = ",".join([str(a.round(3)) for a in com])
        info.append([pdb, compound_name, all_smiles[i], f"pocket_{ith_pocket+1}", com])
info = pd.DataFrame(info, columns=['protein_name', 'compound_name', "smiles", 'pocket_name', 'pocket_com'])

# build the torch dataset:
if os.path.exists(a.dataset_path):
    os.system(f"rm -rf {a.dataset_path}")
    os.system(f"mkdir {a.dataset_path}")
else:
    os.system(f"mkdir {a.dataset_path}")
dataset = TankBind_prediction(a.dataset_path, data=info, protein_dict=protein_dict, compound_dict=compound_dict)
print("Dataset built.")

# load the model:
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("device used: ", device)
logging.basicConfig(level=logging.INFO)
model = get_model(0, logging, device)
model.load_state_dict(torch.load(a.tankbind_model, map_location=device))
model.eval()
model = model.to(device)
dataloader = DataLoader(
    dataset,
    batch_size=64,
    follow_batch=["x","y", "compound_pair"],
    shuffle=False,
    num_workers=32,
)
y_pred_list = []
print(dataset)
print("Running predictions..")
BAs = []
preds = []
for i, data in enumerate(tqdm(dataloader)):
    # try:
        data = data.to(device)
        with torch.no_grad():
            y_pred, affinity_pred = model(data)
        BAs.append(affinity_pred.detach().cpu())
        for i in range(data.y_batch.max() + 1):
            preds.append((y_pred[data['y_batch'] == i]).detach().cpu())

info["affinity"] = torch.cat(BAs)
# save info
info.to_csv(f"{a.output_path}/info_with_BA.csv")
ranked_ba = info.loc[info.groupby(['protein_name', 'smiles'],sort=False)['affinity'].agg('idxmax')].reset_index()
print(ranked_ba)

print("Computing heavy atom coordinates for top 1 BA score for each protein-ligand complex..")
outputs = {}
for i in tqdm(ranked_ba["index"]):
        name = dataset.data.iloc[i].protein_name
        coords = dataset[i].coords
        protein_nodes_xyz = dataset[i].node_xyz
        n_compound = coords.shape[0]
        n_protein = protein_nodes_xyz.shape[0]
        y_pred = preds[i].reshape(n_protein, n_compound)
        y = dataset[i].dis_map.reshape(n_protein, n_compound)
        compound_pair_dis_constraint = torch.cdist(coords, coords)
        pred_dist_info = get_info_pred_distance(coords, y_pred, protein_nodes_xyz,
            compound_pair_dis_constraint,
            n_repeat=1, show_progress=False
        )
        new_coords = pred_dist_info.sort_values("loss")['coords'].iloc[0]
        outputs[name] = new_coords

        # write the ligand's coordinates into a pdb file for inspection:
        savePdb(
            f"{a.nr_lig}/{name.replace('_trimed', '')}.pdb",
            new_coords, 
            f"{a.ligand_output_path}/{name}.pdb"
        )


print(f"Computed {len(outputs)} complexes.")
torch.save(outputs, f"{a.output_path}/outputs.pth")