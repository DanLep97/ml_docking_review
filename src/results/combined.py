from repo import bust
import torch
if __name__ == "__main__":
    import sys
    sys.path.append("../preprocess")
    from qbiolip import load_db
    from repo import get_diffdock_bust_inputs, get_equibind_bust_inputs, get_neuralplexer_bust_inputs
    from repo import get_vina_bust_inputs, get_tankbind_bust_inputs, get_dynamicbind_bust_inputs
    from repo import get_dpl_bust_inputs, get_alphafold_buster_inputs
    from repo import get_crystal_bust_inputs
    from repo import compute_rmsd
    from repo import align_dpl_outputs
    import glob
    import os
    import multiprocessing as mp
    mp.set_start_method("spawn")
    import pandas
    df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
    df = df[df.target]

    # get outputs:
    rec_path = "../../databases/QBioLiP/nr_rec"
    trimmed_rec_path = "../../databases/QBioLiP/nr_rec_trimed"

    dpl_data = "../../databases/DPL/outputs/*"
    dpl_buster_root = "../../databases/DPL/buster_inputs"
    dpl_failed_data = "../../databases/DPL/failed_outputs"
    # align_dpl_outputs(dpl_data, dpl_failed_data) 
    dpl_buster_inputs = get_dpl_bust_inputs(dpl_data, rec_path, trimmed_rec_path, dpl_buster_root)

    diffdock_data = "../../databases/DiffDock/results/user_inference/*"
    diffdock_buster_inputs = get_diffdock_bust_inputs(diffdock_data, rec_path, trimmed_rec_path)

    equibind_data = "../../databases/EquiBind/output/"
    equibind_buster_inputs = get_equibind_bust_inputs(equibind_data, rec_path, trimmed_rec_path)

    neuralplexer_data = "../../databases/NeuralPLexer/inferences/*"
    neuralplexer_buster_inputs = get_neuralplexer_bust_inputs(neuralplexer_data,rec_path, trimmed_rec_path)

    vina_data = "../../databases/Vina/outputs/*"
    vina_buster_inputs = get_vina_bust_inputs(vina_data, rec_path, trimmed_rec_path)

    tankbind_data = "../../databases/TankBind/outputs/outputs.pth"
    tankbind_buster_inputs = get_tankbind_bust_inputs(tankbind_data, rec_path, trimmed_rec_path)

    dynamicbind_data = "/home/daniil/Desktop/Ph.D/DynamicBind/results/*/*"
    dynamicbind_buster_inputs = get_dynamicbind_bust_inputs(dynamicbind_data, rec_path, trimmed_rec_path)
    
    # get the buster and outputs for AF3:
    alphafold_data = "../../databases/alphafold3/outputs/*"
    alphafold_buster_inputs = get_alphafold_buster_inputs(alphafold_data)

    crystal_bust_inputs = get_crystal_bust_inputs(rec_path, df)

# PoseBusters:
    all_buster_inputs = [
        vina_buster_inputs,
        alphafold_buster_inputs,
        equibind_buster_inputs,
        tankbind_buster_inputs,
        diffdock_buster_inputs,
        neuralplexer_buster_inputs,
        dynamicbind_buster_inputs,
        dpl_buster_inputs,
        crystal_bust_inputs
    ]
    with mp.Pool(processes=32) as pool:
        results = []
        for inputs in all_buster_inputs:
            async_results = [
                pool.apply_async(bust, args=item) 
                for item in zip(*inputs)
            ]
            results.append(async_results)
        all_dfs = []
        for async_result_list in results:
            completed_results = []
            for async_result in async_result_list:
                try:
                    # 65 seconds timeout (60 + buffer for file operations)
                    result = async_result.get(timeout=65)
                    if result is not None:
                        completed_results.append(result)
                except mp.TimeoutError:
                    print(f"Task timed out after 65 seconds")
                    continue
                except Exception as e:
                    print(f"Task failed with error: {e}")
                    continue
            if completed_results:
                all_dfs.append(pandas.concat(completed_results))
            else:
                all_dfs.append(pandas.DataFrame())
        vina_bust_df, alphafold_bust_df, equibind_bust_df, tankbind_bust_df, \
        diffdock_bust_df, neuralplexer_bust_df, dynamicbind_bust_df, dpl_bust_df, crystal_df = all_dfs
    alphafold_rmsd = compute_rmsd(alphafold_buster_inputs)
    equibind_rmsd = compute_rmsd(equibind_buster_inputs)
    tankbind_rmsd = compute_rmsd(tankbind_buster_inputs)
    diffdock_rmsd = compute_rmsd(diffdock_buster_inputs)
    neuralplexer_rmsd = compute_rmsd(neuralplexer_buster_inputs)
    dynamicbind_rmsd = compute_rmsd(dynamicbind_buster_inputs)
    dpl_rmsd = compute_rmsd(dpl_buster_inputs)
    vina_rmsd = compute_rmsd(vina_buster_inputs)

    torch.save({
        "bust_df": {
            "AlphaFold3": alphafold_bust_df,
            "EquiBind": equibind_bust_df,
            "TankBind": tankbind_bust_df,
            "DiffDock": diffdock_bust_df,
            "NeuralPLexer": neuralplexer_bust_df,
            "DynamicBind": dynamicbind_bust_df,
            "DPL": dpl_bust_df,
            "Vina": vina_bust_df
        },
        "annex_df": {
            "crystal": crystal_df,
        },
        "rmsds": {
            "AlphaFold3": alphafold_rmsd,
            "EquiBind": equibind_rmsd,
            "TankBind": tankbind_rmsd,
            "DiffDock": diffdock_rmsd,
            "NeuralPLexer": neuralplexer_rmsd,
            "DynamicBind": dynamicbind_rmsd,
            "DPL": dpl_rmsd,
            "Vina": vina_rmsd
        }
    }, "../../databases/results/posebusters_all_df.pt")