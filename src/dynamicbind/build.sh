#!/bin/bash

python prepare_targets.py
cd ../../../DynamicBind/
rm results/* -r
rm data/tmp_* -r
rm data/esm2_output/* -r
rm data/cache_torsion/* -r
DYNAMICBIND_SRC="../ThesisRepo/src/dynamicbind"
declare -a CSV=(`ls $DYNAMICBIND_SRC/csv`)
len=${#CSV[@]}

for ((i=0; i<${len}; i++));
do
	# python run_single_protein_inference.py "$(cat ../ThesisRepo/src/dynamicbind/protein_path.txt)" \
	python run_single_protein_inference.py dummy \
		$DYNAMICBIND_SRC/csv/"${CSV[$i]}" \
	       --savings_per_complex 40 \
	       --protein_path_in_ligandFile \
	       --python /home/daniil/miniconda3/envs/dynamicbind/bin/python \
	       --relax_python /home/daniil/miniconda3/envs/relax/bin/python \
	       --inference_steps 20\
		   --num_workers 20\
	       --header "${CSV[$i]}" 
	mv $DYNAMICBIND_SRC/csv/"${CSV[$i]}" $DYNAMICBIND_SRC/csv_done/"${CSV[$i]}" 
	# rm data/tmp_* -r
	# rm data/esm2_output/* -r
	# rm data/cache_torsion/* -r
done
