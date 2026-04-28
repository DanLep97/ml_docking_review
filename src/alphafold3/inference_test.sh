#!/bin/bash
rm -r ../../databases/alphafold3/multi_chain_outputs/*
rm -r ../../databases/alphafold3/output_test/*

python combine_unique_sequence.py

docker run -it -e JAX_TRACEBACK_FILTERING=off \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/multi_chain_outputs/:/root/af_input/ \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/output_test:/root/af_output \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/databases:/root/public_databases \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/models:/root/models \
    --gpus all \
    alphafold3 \
    python -u run_alphafold.py \
    --input_dir=/root/af_input/ \
    --output_dir=/root/af_output \
    --run_inference=True \