import pandas
import os
import json
import requests
from datetime import datetime
from tqdm import tqdm
import re
import numpy as np

# filter for metals in the binding site:
lig_metal_pattern = r"^(MG|ZN|FE|CU|CO|NI|MN|CA|K|NA|HG|CD|PB|PT|RU|RH|IR|OS|AG|AU)\d{0,1}$"

def load_db(
        qbiolip_db, 
        enzyme_path, 
        build_cols = False,
        no_metal_ligands = False
    ):
    columns = {
        "Q-BioLiP ID": "id",
        "PDB ID": "pdb_id",
        "Assembly ID": "assembly_id",
        "Stoichiometry": "stoichiometry",
        "Resolution (Å)": "resolution",
        "Uniprot ID": "uniprot_id",
        "Assembly Detail":"assembly_detail",
        "Ligand ID":"ligand_id",
        "Site": "site",
        "Ligand Detail": "ligand_detail",
        "Pubmed ID": "pubmed_id",
        "Interaction": "interaction",
        "Relevant":"relevant",
        "Binding Site": "binding_site",
        "Binding Site PDB": "binding_site_pdb",
        "BindingMOAD": "bindingMOAD",
        "PDBbind-CN":"pdbbind_cn",
        "BindingDB": "binding_db"
    }
    if build_cols:
        df = pandas.read_csv(qbiolip_db)
        df.rename(columns=columns, errors="raise", inplace=True)
        # group A: no metal cofactor
        # group B: with metal cofactor outside the binding site
        # group C: with metal cofactor in binding site
        df["group"] = None
        df["lig_canon_smiles"] = None
        df["lig_atoms"] = None
        df["rec_len_resi"] = None
        df["centroid_x"] = None
        df["centroid_y"] = None
        df["centroid_z"] = None
        df["x_size"] = None
        df["y_size"] = None
        df["z_size"] = None
        df["target"] = False
        df["exp_type"] = None
    else:
        df = pandas.read_csv(qbiolip_db, dtype={
            "id": "str",
            "pdb_id": "str",
            "assembly_id": "str",
            "stoichiometry": "str",
            "resolution": np.float32,
            "uniprot_id": "str",
            "assembly_detail": "str",
            "ligand_id": "str",
            "site": "str",
            "ligand_detail": "str",
            "pubmed_id": "str",
            "interaction": "str",
            "relevant": "str",
            "binding_site": "str",
            "binding_site_pdb": "str",
            "bindingMOAD": "str",
            "pdbbind_cn": "str",
            "binding_db": "str",
            "group": "str",
            "lig_canon_smiles": "str",
            "lig_atoms": "str",
            "rec_len_resi": np.float32,
            "centroid_x": np.float32,
            "centroid_y": np.float32,
            "centroid_z": np.float32,
            "x_size": np.float32,
            "y_size": np.float32,
            "z_size": np.float32,
            "target": np.bool_,
            "exp_type": "str",
            "size_x": np.float32,
            "size_y": np.float32,
            "size_z": np.float32
        }, low_memory=False)
    # pdb with annotated catalytic residues are assumed to be enzymes
    enzyme_file = open(enzyme_path, "r")
    enzyme_uniprot_ids = [uniprot_id.replace("\n", "") for uniprot_id in enzyme_file]
    # retreive cofactors:
    base_cofactors = [
        "TDP", "FAD", "FMN", "NAD", "PNS", "COA"
        "PLP", "GSH", "BTN", "FFO", "B12", "ASC",
        "MQ7", "UQ1", "MGD", "H4B", "MDO", "SAM",
        "F43", "COM", "TP7", "HEA", "DPM", "PQQ",
        "TPQ", "TRQ", "LPA", "GOL", "ATP", "ADP",
        "UNK"
    ]
    cofactors = [v[0]["cofactors"] for _,v in json.load(open("../../databases/cofactors.json")).items()]
    cofactors = [c for c_l in cofactors for c in c_l if c is not None] + base_cofactors
    # retreive crystal artefacts:
    xtal_artefacts = [l.replace("\n", "") for l in open("../../databases/crystal_artefacts.txt")]
    xtal_artefacts = [l for l in xtal_artefacts if re.match(lig_metal_pattern, l) == None]
    df = df.loc[
        (df.interaction == "yes") & (df.relevant == "yes")
        & (df.ligand_id != "DNA") & (df.ligand_id != "III") #DNA and peptides are out
        & (df.ligand_id != "DNA+RNA") # DNA interacting with RNA are out too
        & (df.ligand_id != "RNA") # remove RNA ligands
        & (df.uniprot_id.isin(enzyme_uniprot_ids)) # focus only on enzymes
        & (df.ligand_id != "") # ligand ID should be provided
        & (~df.ligand_id.isin(cofactors)) # focus only on true ligands
        & (~df.ligand_id.isna()) # same as earlier, ligand ID should be provided
        & (~df.ligand_id.isin(xtal_artefacts))
    ]
    if not build_cols:
        df = df.loc[~df.lig_canon_smiles.isna()]
    if no_metal_ligands:
        df = df.loc[~df.ligand_id.str.match(lig_metal_pattern)]
    print(f"Len of the database: {df.shape[0]}")
    print(f"Columns: {', '.join(df.columns)}")
    df.reset_index(inplace=True, drop=True)
    return df # only protein-small ligands and protein-metal interaction are kept

aa_3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    "SEC": "U"
}

aa_1to3 = {value:key for key,value in aa_3to1.items()}

if __name__ == "__main__":
    import argparse
    import torch

    arg_parser = argparse.ArgumentParser(description="""
        Reads the protein-ligand annotation file from Q-BioLiP database.
        Website where the .list file can be downloaded: https://yanglab.qd.sdu.edu.cn/Q-BioLiP/Download/
        The file is under "interaction-based data" -> Protein Ligand -> annotation (non redudant).
        The file extension is .list which should be modified to .csv.
        The "load_db" function returns the pandas DataFrame object.
        This scripts annotates each receptor into 3 categories: with metal ions in BS, no metal ions, 
        and metal ions outside the protein
    """)

    arg_parser.add_argument(
        "--cid",
        help="""
        Path to the cid.csv file containing annotations between cid and SMILES.
        Default: ../../databases/components.csv 
        """,
        default="../../databases/components.csv"

    )
    arg_parser.add_argument(
        "--ec-uniprot-ids",
        help="""
        Path to the txt file containing the list of enzyme uniprot ids.
        This text file is generated using `src/preprocess/brenda_json_ec.py file.
        Default: ../../databases/ec_uniprot_ids.txt 
        """,
        default="../../databases/ec_uniprot_ids.txt"
    )
    arg_parser.add_argument(
        "--csv-path",
        help="""
        Path to the annotation file. Default: ../../databases/QBioLiP/qbiolip.csv
        """,
        default="../../databases/QBioLiP/qbiolip.csv" 
    )
    arg_parser.add_argument(
        "--receptor-path",
        help="""
        Path to receptor folder. Default: ../../databases/QBioLiP/nr_rec
        """,
        default="../../databases/QBioLiP/nr_rec"
    )
    arg_parser.add_argument(
        "--ligand-path",
        help="""
        Path to receptor folder. Default: ../../databases/QBioLiP/nr_lig
        """,
        default="../../databases/QBioLiP/nr_lig"
    )

    a = arg_parser.parse_args()

    df = load_db(a.csv_path, a.ec_uniprot_ids, build_cols=True)
    failed_cases = 0

    # get the smiles of CID ids using cid.csv file
    components = pandas.read_csv(a.cid)
    cids = df.ligand_id.unique()
    three_letter_codes = components.three_letter_code.tolist()
    print("Annotating ligands with canonical SMILES representations..")
    for cid in tqdm(cids):
        if cid in three_letter_codes:
            canonical_smiles = components.loc[components.three_letter_code == cid].canonical_smiles.values[0]
            if isinstance(canonical_smiles, str):
                df.loc[df.ligand_id == cid, ["lig_canon_smiles"]] = canonical_smiles
    annotated_cases = df.loc[df.lig_canon_smiles != "unknown"].ligand_id.unique().shape[0]
    perc = annotated_cases/df.ligand_id.unique().shape[0]*100
    print(f"{perc:.2f}% of unique ligand ids were annotated with SMILES representations")

    # compute number of atoms in ligands and number of residues in assembly:
    def num_atoms(x):
        if x is None:
            return None
        return "".join(re.findall(r"[A-Z]*", x))
    def num_res(assembly):
        rec_path = f"{a.receptor_path}/{assembly}.pdb"
        if os.path.exists(rec_path):
            with open(rec_path, "r") as rec_f:
                resi = set(int(l[22:26].strip()) for l in rec_f if l[0:6].strip() == "ATOM")
                return len(resi)
        else:
            return None
    print("Reading number of atoms per ligands..")
    df["lig_atoms"] = df.lig_canon_smiles.apply(lambda x: num_atoms(x))
    print("Reading number of receptor residues..")
    df["rec_len_resi"] = df.assembly_id.apply(lambda x: num_res(x))

    # get the experiment type for each pdb:
    def get_experiment_type(pdb_id):
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            method = data.get('rcsb_entry_info', {}).get('experimental_method', 'Unknown')
            return method
        else:
            return None
    # print("Fetching experiment type for each pdb_ids:")
    # pdb_ids = df.pdb_id.unique()
    # for pdb_id in tqdm(pdb_ids):
    #     exp_type = get_experiment_type(pdb_id)
    #     df.loc[df.pdb_id == pdb_id, ["exp_type"]] = exp_type
    # print("Experimental techniques used to generate structure:", df.exp_type.unique())
    

    for i in range(df.shape[0]):
        try:
            case = df.iloc[i]
            # if case.exp_type == "NMR" or case.exp_type is None:
            #     continue
            if case.lig_canon_smiles is None:
                continue
            lig_atoms = "".join(re.findall(r"[A-Za-z]*", case.lig_canon_smiles))
            # only protein-ligand interactions and ligands with more than 5 atoms are of interest:
            if re.match(lig_metal_pattern, case.ligand_id) or len(lig_atoms) < 5:
                continue
            ligand_file = f"{a.ligand_path}/{case.ligand_detail}.pdb"
            if not os.path.exists(ligand_file):
                continue
            with open(ligand_file) as lig_f:
                elements = set(l[76:78].strip() for l in lig_f if l[0:6].strip() == "HETATM")
                if any(re.match(lig_metal_pattern, el) for el in elements):
                    continue

            receptor_file = f"{a.receptor_path}/{case.assembly_id}.pdb"
            assert os.path.isfile(receptor_file), f"Receptor file {receptor_file} does not exist."

            br_strs = case.binding_site_pdb.split("  ") # separate the binding residues string per chain
            assert len(br_strs) > 0, f"No binding residues for {receptor_file} with {case.ligand_id}"

            chains = [c.upper() for c in re.findall(r"[A-Za-z](?=:)", "  ".join(br_strs))] # get the chains involved in BS
            res = [[aa_1to3[r.upper()] for r in re.findall(r"[A-Za-z](?=[0-9]{1,})", s)] for s in br_strs] # get the residues involved in BS
            resi = [[r for r in re.findall(r"\d+", s)] for s in br_strs] # get the residue number involved in BS
            ca_info = [
                (chain, res[i][j], resi[i][j]) 
                for i, chain in enumerate(chains)
                for j in range(len(res[i])) 
            ]

            # open the pdb file of the receptor, read only carbon alpha ATOM records:
            atom_rec_records = [l.replace("\n", "") for l in open(receptor_file, "r") if\
                l[0:6].strip() == "ATOM" 
            ]
            atom_lig_records = [l.replace("\n", "") for l in open(ligand_file, "r") if\
                l[0:6].strip() == "HETATM"
            ]
            assert len(atom_rec_records)>0, f"No ATOM records available in {receptor_file}"
            assert len(atom_lig_records)>0, f"No HETATM records available in {ligand_file}"
            # get alpha carbon atom records for binding residues, aggregated per chain:
            ca = [l for chain,res,resi in ca_info for l in atom_rec_records if l[21] == chain
                    and l[22:26].strip() == resi
                    and l[12:16].strip() == "CA"
                    and l[17:20].strip() == res]

            # verify all alpha carbon are retreived:
            assert len(ca_info) == len(ca),\
            f"Number of CA atoms discrepancy in {receptor_file}: expected {len(ca_info)}, {len(ca)} gathered."

            # get xyz coordinates of the binding site receptor alpha carbon atoms:
            xyz_rec = [
                [
                    float(r[30:38].strip()), #x
                    float(r[38:46].strip()), #y
                    float(r[46:54].strip())  #z
                ] for r in ca
            ]
            xyz_lig = [
                [
                    float(r[30:38].strip()), #x
                    float(r[38:46].strip()), #y
                    float(r[46:54].strip())  #z
                ] for r in atom_lig_records
            ]
            xyz_rec = torch.tensor(xyz_rec) # convert to tensor for geometric manipulation
            xyz_lig = torch.tensor(xyz_lig)
            # print(xyz)
            assert xyz_rec.shape[0] != 0, f"No xyz coordinates available for carbon alpha ATOM record from {receptor_file}"
            assert xyz_lig.shape[0] != 0, f"No xyz coordinates available for HETATM record from {ligand_file}"

            # geometrically define the binding site:
            # the centroid is the center of the receptor binding site:
            centroid = xyz_rec.mean(dim=0) # the center of the binding site is the mean value of CA xyz coordinates

            # the size is defined by the max distance between ligand atoms, on each axis:
            # from vina paper, it is recommended that the size for each axis is size of ligand (on the axis) + 15A
            # or higher than 22.5 total.
            xyz_lig_expd = xyz_lig.expand(xyz_lig.shape[0], xyz_lig.shape[0], -1)
            axis_dist = (xyz_lig_expd - xyz_lig_expd.transpose(0,1)).abs().max(dim=0).values.max(dim=0).values + 15
            axis_dist[axis_dist<22.5] = 22.5 # minimum value of 22.5

            # define the binding site grid for A,B,C categorization:
            xyz_rec_expd = xyz_rec.expand(xyz_rec.shape[0], xyz_rec.shape[0], -1)
            rec_axis_dist = (xyz_rec_expd - xyz_rec_expd.transpose(0,1)).abs().max(dim=0).values.max(dim=0).values + 10
            grid_boundaries = torch.stack([centroid - 1/2*rec_axis_dist, centroid + 1/2*rec_axis_dist])

            # gather all metal ligand files from this receptor except the current ligand:
            met_cof = df.loc[
                (df.ligand_detail.str.contains(case.assembly_id))
                & (df.ligand_id.str.match(lig_metal_pattern))
            ]
            met_cof_pdb = [f"{a.ligand_path}/{l}.pdb" for l in met_cof.ligand_detail]
            met_cof_pdb_exist = all(os.path.isfile(f) for f in met_cof_pdb)
            assert met_cof_pdb_exist, f"Not all metal ion pdb file are available for {receptor_file}"
            # get the list of all hetatm records for each chain, used to filter for metals:
            metal_coords = []
            for met_pdb in met_cof_pdb:
                metal_coord = [(
                    float(l[30:38].strip()), #x
                    float(l[38:46].strip()), #y
                    float(l[46:54].strip())  #z
                ) for l in open(met_pdb, "r") if l[:6].strip() == "HETATM"]
                assert len(metal_coord)==1, f"Metal ion file {met_pdb} has more than one hetatm record"
                metal_coords.extend(metal_coord)
            metal_coords = torch.tensor(metal_coords)
            assert metal_coords.shape[0] == len(met_cof_pdb), f"Not all metal coordinates are available for {receptor_file}"
            # add geometric binding site information:
            df.loc[i, [
                "centroid_x", "centroid_y", "centroid_z",
                "size_x", "size_y", "size_z"
            ]] = [
                centroid[0].item(), centroid[1].item(), centroid[2].item(),
                axis_dist[0].item(), axis_dist[1].item(), axis_dist[2].item()
            ]
            
            if metal_coords.shape[0] > 0:
                in_bounds = (
                    (metal_coords > grid_boundaries[0]) & (metal_coords < grid_boundaries[1])
                ).nonzero(as_tuple=True)
                # list of all metal atoms indices which are in bounds:
                in_bounds = [i for i in range(metal_coords.shape[0]) if (i == in_bounds[0]).nonzero().shape[0] == 3]
                if len(in_bounds) > 0:
                    df.loc[i, "group"] = "C"
                else:
                    df.loc[i, "group"] = "B"
            else:
                df.loc[i, "group"] = "A"

        except AssertionError as e:
            print(e)
            failed_cases=failed_cases+1
        except TypeError as e:
            print(e)
            failed_cases=failed_cases+1
        except ValueError as e:
            print(e)
            failed_cases = failed_cases+1

    # print some stats:
    print("Failed cases:",failed_cases)
    print(f"Percentage: {failed_cases/df.shape[0]*100:.3f}")

    # save csv:
    df.to_csv(a.csv_path, index=False)
    print(f"csv saved at {a.csv_path} with {df.shape[0]} rows")