from rdkit import Chem
from rdkit.Chem import AllChem
from qbiolip import load_db

df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
df = df[df.target]

n_non_unique = 0
names = []
for i in range(len(df)):
    row = df.iloc[i]
    src_f = f"/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_lig/{row.ligand_detail}_randconf.pdb" 
    with open(src_f, "r") as f:
        atoms = [l for l in f if l[:6].strip() == "HETATM"]
    atom_names = [a[12:16].strip() for a in atoms]
    atom_names_unique = set(atom_names)
    if len(atom_names) != len(atom_names_unique):
        n_non_unique+=1
        names.append(row.ligand_detail)
print("Number of n_non_unique:", n_non_unique)
print("Names:", names)
