#!/bin/bash

# don't forget to activate conda env before running this script!
python prepare_targets.py
# rm -rf ../../databases/NeuralPLexer/inferences/*

declare -a RECEPTORS=(`cat rec_paths.txt`)
declare -a LIGANDS=(`cat lig_paths.txt`)
declare -a OUTPUTS=(`cat output_paths.txt`)

len=${#OUTPUTS[@]}
for (( i=0; i<${len}; i++));
do
	printf "Current index %d with value %s\n" $i "${OUTPUTS[$i]}"
	mkdir "${OUTPUTS[$i]}"

	neuralplexer-inference --task=batched_structure_sampling \
		--input-receptor "${RECEPTORS[$i]}" \
		--input-ligand "${LIGANDS[$i]}" \
		--use-template \
		--input-template "${RECEPTORS[$i]}" \
		--out-path "${OUTPUTS[$i]}" \
		--model-checkpoint /home/daniil/Desktop/Ph.D/NeuralPLexer/models/complex_structure_prediction.ckpt \
		--n-samples 10 \
		--chunk-size 1 \
		--num-steps=40 \
		--cuda \
		--sampler=langevin_simulated_annealing\
		--rank-outputs-by-confidence
done
