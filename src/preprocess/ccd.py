import argparse
import pandas
from pdbecif.mmcif_io import CifFileReader

arg_parser = argparse.ArgumentParser(
    description="""
    Parse a mmCIF file into a csv containing the type and three letter for each entry.
    The csv file is saved in the same location as the mmCIF file.
    """
)

arg_parser.add_argument(
    "--cif",
    help="Path to mmCIF file. Default: ../../databases/components.cif",
    default="../../databases/components.cif"
)

a = arg_parser.parse_args()
category_names = ["_chem_comp", "_pdbx_chem_comp_descriptor"]
db = CifFileReader().read(a.cif, only=category_names)

ids = []
names = []
three_letter_codes = []
canonical_smiles = []
types = []

for data_id, descr in db.items():
    # Check that each data points have all categories
    id, name, three_letter_code, smiles, t = "unknown", None, None, None, None
    if "_chem_comp" in descr.keys():
        id = descr["_chem_comp"]["id"]
        name = descr["_chem_comp"]["name"]
        three_letter_code = descr["_chem_comp"]["three_letter_code"]
        t = descr["_chem_comp"]["type"]
    if "_pdbx_chem_comp_descriptor" in descr.keys():
        # fetch canonical smiles representation
        smiles_types = descr["_pdbx_chem_comp_descriptor"]["type"]
        if "SMILES_CANONICAL" in smiles_types:
            canonical_idx = smiles_types.index("SMILES_CANONICAL")
            smiles = descr["_pdbx_chem_comp_descriptor"]["descriptor"][canonical_idx]
    ids.append(id)
    names.append(name)
    three_letter_codes.append(three_letter_code)
    canonical_smiles.append(smiles)
    types.append(t)
    
        
df = pandas.DataFrame(
    {
        "id":ids,
        "name":names,
        "three_letter_code":three_letter_codes,
        "canonical_smiles":canonical_smiles, # canonical only!
        "type": types
    }
)

csv_path = a.cif.replace(".cif", ".csv")
df.to_csv(csv_path, index=False)
print(f"csv file saved at {csv_path}")
