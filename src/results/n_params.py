import torch

equibind_model_p = "/home/daniil/Desktop/Ph.D/EquiBind/runs/flexible_self_docking/best_checkpoint.pt" # equibind
tankbind_model_p = "/home/daniil/Desktop/Ph.D/TankBind/saved_models/self_dock.pt"
diffdock_model_p = [
    "/home/daniil/Desktop/Ph.D/DiffDock/workdir/v1.1/confidence_model/best_model_epoch75.pt",
    "/home/daniil/Desktop/Ph.D/DiffDock/workdir/v1.1/score_model/best_ema_inference_epoch_model.pt"
]
neuralplexer_model_p = "/home/daniil/Desktop/Ph.D/NeuralPLexer/models/complex_structure_prediction.ckpt"
dynamicbind_model_p = "/home/daniil/Desktop/Ph.D/DynamicBind/workdir/big_score_model_sanyueqi_with_time/pro_ema_inference_epoch138_model.pt"
dpl_model_p = "/home/daniil/Desktop/Ph.D/DiffusionProteinLigand/checkpoints/DPLS_v1.ckpt"

n_params = []

equibind_model = torch.load(equibind_model_p)["model_state_dict"]
n_params.append({
    "name": "EquiBind",
    "n":sum([p.numel() for p in equibind_model.values()])
})

tankbind_model = torch.load(tankbind_model_p)
n_params.append({
    "name": "TankBind",
    "n": sum([p.numel() for p in tankbind_model.values()])
})

diffdock_n_params = 0
for m in diffdock_model_p:
    model = torch.load(m)
    diffdock_n_params += sum([p.numel() for p in model.values()])
n_params.append({
    "name": "DiffDock",
    "n": diffdock_n_params
})

neuralplexer_model = torch.load(neuralplexer_model_p)["state_dict"]
n_params.append({
    "name": "NeuralPLexer",
    "n": sum([p.numel() for p in neuralplexer_model.values()])
})

dynamicbind_model = torch.load(dynamicbind_model_p)
n_params.append({
    "name": "DynamicBind",
    "n":sum([p.numel() for p in dynamicbind_model.values()])
})

dpl_model = torch.load(dpl_model_p)["state_dict"]
n_params.append({
    "name": "DiffusionProteinLigand",
    "n":sum([p.numel() for p in dpl_model.values()])
})

for m in n_params:
    print(f"For model {m['name']}, number of parameters: {m['n']:,}")