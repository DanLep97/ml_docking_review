#!/bin/bash
declare -a JSON_PATHS=(`ls /root/af_input/*.json`)
declare -a JSONS=()
for path in "${JSON_PATHS[@]}"; do
    # Extract filename using basename
    filename=$(basename "$path")
    # Add filename to the array
    JSONS+=("$filename")
done
len=${#JSONS[@]}
for ((i=0; i<${len}; i++));
do
        python -u run_alphafold.py \
        --json_path=/root/af_input/"${JSONS[$i]}" \
        --output_dir=/root/af_output \
        --run_inference=True \
        --jax_compilation_cache_dir . \
        --norun_data_pipeline
done