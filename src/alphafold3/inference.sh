#!/bin/bash
rm -r ../../databases/alphafold3/multi_chain_outputs/*
#rm -r ../../databases/alphafold3/outputs/*
rm ../../databases/alphafold3/cif/*
python combine_unique_sequence.py

# docker run -it -e JAX_TRACEBACK_FILTERING=off \
docker run -it --privileged \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/multi_chain_outputs/:/root/af_input/ \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/outputs:/root/af_output \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/databases:/root/public_databases \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/databases/alphafold3/models:/root/models \
    --volume /home/daniil/Desktop/Ph.D/ThesisRepo/src/alphafold3/docker:/root/custom_script/ \
    --gpus all \
    alphafold3 \
    /root/custom_script/run_alphafold.sh