from qbiolip import load_db
import glob
from tqdm import tqdm
import re
import torch
from pprint import pprint

rec_path = "../../databases/QBioLiP/nr_rec"
lig_path = "../../databases/QBioLiP/nr_lig"
# pattern_text = "MG|ZN|FE|CU|CO|NI|MN|CA|K|NA|HG|CD|PB|PT|RU|RH|IR|OS|AG|AU"
# lig_metal_pattern = r"^(MG|ZN|FE|CU|CO|NI|MN|CA|K|NA|HG|CD|PB|PT|RU|RH|IR|OS|AG|AU)\d{0,1}$"
pattern_text = "PD"
lig_metal_pattern = r"^(PD)\d{0,1}$"
metals = pattern_text.split("|")

# Create dictionary
metal_dict = {metal: 0 for metal in metals}

interaction_threshold = 4
all_ligand_paths = glob.glob(f"{lig_path}/*.pdb")
all_ligand_CCDs = [p.split("/")[-1].split("_")[-2] for p in all_ligand_paths]
metal_ligand_paths = [
    all_ligand_paths[i] for i in range(len(all_ligand_paths))
    if re.search(lig_metal_pattern, all_ligand_CCDs[i])
]

for i,ligand_file in enumerate(metal_ligand_paths):
    metal_ligand_filename = ligand_file.split("/")[-1]
    print(metal_ligand_filename)
    metal_type = all_ligand_CCDs[i][:2]
    if "8g27" in metal_ligand_filename or metal_type not in list(metal_dict.keys()):
        continue
    assembly_id = "_".join(metal_ligand_filename.split("_")[:2])
    receptor_file = f"{rec_path}/{assembly_id}.pdb"

    atom_rec_records = [l.replace("\n", "") for l in open(receptor_file, "r") if\
        l[0:6].strip() == "ATOM" 
    ]
    atom_lig_records = [l.replace("\n", "") for l in open(ligand_file, "r") if\
        l[0:6].strip() == "HETATM"
    ]
    xyz_rec = torch.tensor([
        [
            float(r[30:38].strip()), #x
            float(r[38:46].strip()), #y
            float(r[46:54].strip())  #z
        ] for r in atom_rec_records
    ])
    xyz_lig = torch.tensor([
        [
            float(r[30:38].strip()), #x
            float(r[38:46].strip()), #y
            float(r[46:54].strip())  #z
        ] for r in atom_lig_records
    ])
    if xyz_lig.shape[0]!=1 or metal_type not in list(metal_dict.keys()) or metal_ligand_filename == "8g27_1_MG_N.pdb":
        continue
    distogram = (xyz_rec - xyz_lig).pow(2).sum(dim=1).sqrt()
    in_threshold = (distogram < interaction_threshold).nonzero().shape[0]>0
    if in_threshold:
        metal_dict[metal_type] = metal_dict[metal_type] + 1
    
pprint(metal_dict)