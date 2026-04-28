#!/bin/bash
rm -r ../../databases/alphafold3/inputs/*
rm -r ../../databases/alphafold3/unique_chain_outputs/*
python prepare_targets.py

# Run this script first. Then, run inference.sh
# This script is used to compute the MSA for every unique protein chain.
docker run -it\
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/inputs:/root/af_input/ \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/unique_chain_outputs:/root/af_output \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/databases:/root/public_databases \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/models:/root/models \
    --gpus all \
    alphafold3 \
    python -u run_alphafold.py \
    --input_dir=/root/af_input/ \
    --output_dir=/root/af_output \
    --run_inference=False \
