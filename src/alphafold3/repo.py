from rdkit import Chem
from rdkit.Chem import AllChem
import os

def pdb_to_ccd_mmcif(pdb_file, output_cif, comp_id="UNL"):
    """
    Convert PDB to CCD-like mmCIF using RDKit for accurate bond/charge detection.
    
    Args:
        pdb_file (str): Input PDB file path.
        output_cif (str): Output mmCIF file path.
        comp_id (str): Chemical component ID (default: "UNL").
    """
    # Load PDB with RDKit
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("RDKit failed to parse PDB file")

    # Sanitize and assign properties
    Chem.SanitizeMol(mol)
    AllChem.ComputeGasteigerCharges(mol)  # For partial charges (optional)
    
    # ===== (1) Singular Fields =====
    cif_data = {
        "_chem_comp": {
            "id": comp_id,
            "name": "?",
            "type": "non-polymer",
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "mon_nstd_parent_comp_id": "?",
            "pdbx_synonyms": "?",
            "formula_weight": f"{Chem.rdMolDescriptors.CalcExactMolWt(mol):.2f}",
        },
        "_chem_comp_atom": [],
        "_chem_comp_bond": [],
    }

    # ===== (2) Per-Atom Fields =====
    for atom in mol.GetAtoms():
        cif_data["_chem_comp_atom"].append({
            "comp_id": comp_id,
            "atom_id": atom.GetPDBResidueInfo().GetName().strip(),
            "type_symbol": atom.GetSymbol(),
            "charge": atom.GetFormalCharge(),
            "pdbx_leaving_atom_flag": "N",
            "pdbx_model_Cartn_x_ideal": mol.GetConformer().GetAtomPosition(atom.GetIdx()).x,
            "pdbx_model_Cartn_y_ideal": mol.GetConformer().GetAtomPosition(atom.GetIdx()).y,
            "pdbx_model_Cartn_z_ideal": mol.GetConformer().GetAtomPosition(atom.GetIdx()).z,
        })

    # ===== (3) Per-Bond Fields =====
    for bond in mol.GetBonds():
        cif_data["_chem_comp_bond"].append({
            "atom_id_1": mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetPDBResidueInfo().GetName().strip(),
            "atom_id_2": mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetPDBResidueInfo().GetName().strip(),
            "value_order": {
                1: "SING",
                2: "DOUB",
                3: "TRIP",
                1.5: "AROM",
            }.get(bond.GetBondTypeAsDouble(), "SING"),
            "pdbx_aromatic_flag": "Y" if bond.GetIsAromatic() else "N",
        })

    # ===== (4) Write mmCIF File =====
    with open(output_cif, "w") as f:
        # Write header
        f.write(f"data_{comp_id}\n#\n")
        
        # Write _chem_comp block
        for k, v in cif_data["_chem_comp"].items():
            f.write(f"_chem_comp.{k}\t{v}\n")
        f.write("#\n")
        
        # Write _chem_comp_atom block
        f.write("loop_\n")
        f.write("\n".join(f"_chem_comp_atom.{k}" for k in cif_data["_chem_comp_atom"][0].keys()) + "\n")
        for atom in cif_data["_chem_comp_atom"]:
            f.write("\t".join(str(v) for v in atom.values()) + "\n")
        f.write("#\n")
        
        # Write _chem_comp_bond block
        f.write("loop_\n")
        f.write("\n".join(f"_chem_comp_bond.{k}" for k in cif_data["_chem_comp_bond"][0].keys()) + "\n")
        for bond in cif_data["_chem_comp_bond"]:
            f.write("\t".join(str(v) for v in bond.values()) + "\n")
        f.write("#\n")

if __name__ == "__main__":
    # Example usage
    pdb_to_ccd_mmcif(
        "/home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_lig/1gui_1_BGC_B_randconf.pdb", 
        "./output.cif", 
        comp_id="LIG"
    )