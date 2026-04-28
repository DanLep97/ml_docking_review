import torch
import pandas
import random
import re
import sys
sys.path.append("/home/daniil/Desktop/Ph.D/ThesisRepo/src/preprocess")
from qbiolip import load_db
import re

df = load_db("../../databases/QBioLiP/qbiolip.csv", "../../databases/ec_uniprot_ids.txt")
state = df.loc[df.ligand_detail == "1kdg_2_MAN_Q1"].iloc[0]
pdb = f"../../databases/QBioLiP/nr_rec/{state.assembly_id}.pdb"
print("Working with: ", pdb)
centroid = torch.tensor([
    state.centroid_x, state.centroid_y, state.centroid_z
])
size = torch.tensor([
    state.size_x, state.size_y, state.size_z
])
print("centroid:", centroid)
print("size:", size)

sphere_dist = 2
sphere_radius = .5

depth = (size/sphere_dist).ceil().int()
x = torch.tensor([centroid[0]-1/2*size[0] + i*sphere_dist for i in range(depth[0])])
y = torch.tensor([centroid[1]-1/2*size[1] + i*sphere_dist for i in range(depth[1])])
z = torch.tensor([centroid[2]-1/2*size[2] + i*sphere_dist for i in range(depth[2])])
grid_x,grid_y,grid_z = torch.meshgrid(x,y,z,indexing="xy")
grid = torch.stack([grid_x, grid_y, grid_z], dim=3).reshape(-1,3)
print(f"Grid shape: {grid_x.shape}, number of spheres: {grid.shape[0]}")
print(f"Pdb used: {pdb}")
sphere_list = [
    6.0, 1.0, 0, 0,
    7.0, *centroid.tolist(), 1 #start with sphere in the center
]
for xyz in grid:
    xyz=xyz.tolist()
    sphere_list.extend([6.0, 0.8, 0, 0])
    sphere_list.extend([7.0, *xyz, sphere_radius])

template = f"""
from pymol import cmd
from pymol.cgo import *

pdb = "{pdb}"
cmd.load(pdb)

sphere_list = {sphere_list}
# for i,s in enumerate(sphere_list):
#     cmd.load_cgo(s, "sphere_"+str(i))
cmd.load_cgo(sphere_list, "spheres", 1)
"""

f = open("grid.py", "w")
f.write(template)
f.close()
