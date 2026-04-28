#!/bin/bash

# don't forget to activate conda env before running this script!
python prepare_targets.py
# rm -rf ../../databases/DPL/*/*

declare -a RECEPTORS=(`cat rec_paths.txt`)
declare -a LIGANDS=(`cat lig_paths.txt`)
declare -a OUTPUTS=(`cat output_paths.txt`)
DPL=/home/daniil/Desktop/Ph.D/DiffusionProteinLigand
CKPT=/home/daniil/Desktop/Ph.D/DiffusionProteinLigand/checkpoints/DPLS_v1.cpkt

len=${#OUTPUTS[@]}
for (( i=0; i<${len}; i++));
do
	printf "Current index %d with value %s\n" $i "${OUTPUTS[$i]}"

	python "${DPL}/generate.py" \
		--ckpt_path "${DPL}/checkpoints/DPLS_v1.ckpt" \
		--output_dir "${OUTPUTS[$i]}" \
		--protein "${RECEPTORS[$i]}" \
		--ligand "${LIGANDS[$i]}" \
		--num_samples 8 \
		--num_steps 64
done