import os
import shutil
from glob import glob

if __name__=="__main__":
    import sys
    sys.path.append("../preprocess")
    import argparse
    from qbiolip import load_db
    import multiprocessing as mp
    import re
    arg_parser = argparse.ArgumentParser(description="""
        Using SMILES representations from the annotation file, generates a pdbqt string file and a receptor file
        for each case.
    """)
    arg_parser.add_argument(
        "--ntargets",
        help="Number of targets to randomly select from each class. Default: 30",
        default=30,
        type=int
    )
    arg_parser.add_argument(
        "--qbiolip",
        help="Path to the qbiolip annotation csv file. Default: ../../databases/QBioLiP/qbiolip.csv",
        default="../../databases/QBioLiP/qbiolip.csv"
    )
    arg_parser.add_argument(
        "--ec-uniprot-ids",
        help="Path to the biolip database. Used to filter enzymes. Default: ../../databases/BioLiP/BioLiP_nr.csv",
        default="../../databases/ec_uniprot_ids.txt"
    )
    arg_parser.add_argument(
        "--data",
        help="""
        Path to EquiBind data folder to generate required file structure.
    """,
        default="../../databases/EquiBind/data"
    )
    arg_parser.add_argument(
        "--lig-out",
        help="""
        Path to EquiBind output containing the .sdf files.
    """,
        default="../../../EquiBind/data/results/output"
    )
    arg_parser.add_argument(
        "--lig",
        help="""
        Path to qbiolip ligand files. Default: ../../databases/QBioLiP/nr_lig
    """,
        default="../../databases/QBioLiP/nr_lig"
    )
    arg_parser.add_argument(
        "--rec",
        help="""
        Path to qbiolip receptor files. Default: ../../databases/QBioLiP/nr_rec
    """,
        default="../../databases/QBioLiP/nr_rec"
    )

    a = arg_parser.parse_args()
    df = load_db(a.qbiolip, a.ec_uniprot_ids)
    df = df.loc[df.target]
    df.reset_index(inplace=True, drop=True)
    modelled_paths = glob(f"{a.lig_out}/*")
    modelled = [p.split("/")[-1] for p in modelled_paths]
    print(modelled[:10])
    print("Number of complexes:", df.shape)
    for ind in df.index:
        case = df.loc[ind]
        lig_detail = case.ligand_detail
        data_folder = f"{a.data}/{lig_detail}"
        data_folder_trimed = f"{data_folder}_trimed"
        try:
            if case.ligand_detail not in modelled:
                os.mkdir(data_folder)
            if case.ligand_detail + "_trimed" not in modelled:
                os.mkdir(data_folder_trimed)
        except FileExistsError:
            continue
        ligand_src = f"{a.lig}/{lig_detail}_randconf.pdb"
        receptor_src = f"{a.rec}/{case.assembly_id}.pdb"
        receptor_trimed_src = f"{a.rec}_trimed/{case.ligand_detail}.pdb"

        ligand_dest = f"{data_folder}/{lig_detail}_ligand.pdb"
        ligand_trimed_dest = f"{data_folder_trimed}/{lig_detail}_trimed_ligand.pdb"
        receptor_dest = f"{data_folder}/{lig_detail}_protein.pdb"
        receptor_trimed_dest = f"{data_folder_trimed}/{lig_detail}_trimed_protein.pdb"

        if case.ligand_detail not in modelled:
            shutil.copyfile(receptor_src, receptor_dest)
            shutil.copyfile(ligand_src, ligand_dest)
        if case.ligand_detail + "_trimed" not in modelled:
            shutil.copyfile(receptor_trimed_src, receptor_trimed_dest)
            shutil.copyfile(ligand_src, ligand_trimed_dest)
    ncases = len(glob(f"{a.data}/*"))
    print(f"Number of targets in the database: {df.shape[0]}")
    print(f"Number of prepared targets in {a.data}: {ncases}")