#!/bin/bash

# use this script to generate docked structures using EquiBind github repository
# available at https://github.com/HannesStark/EquiBind/tree/main
# To run this script, first modify the inference config file as explained in 
# the readme file of the repository. Before using this script, the python script 
# /src/equibind/prepare_targets.py has to build the correct folder structure
# required to run this bash file.
# The conda environment required for EquiBind has to be activated before running
# this script.

# Path to the repository
#EQUIBIND_PATH=/home/daniil/Desktop/Ph.D/EquiBind
EQUIBIND_PATH=/home/daniil/Downloads/EquiBind
# Path to input data
INPUT_PATH=/home/daniil/Desktop/Ph.D/ThesisRepo/databases/EquiBind/data

# data preprocessing:
rm $INPUT_PATH/* -rf
python prepare_targets.py
cd $EQUIBIND_PATH
# python ./data_preparation/openbabel_receptors.py
# python ./data_preparation/reduce_receptors.py
# rm $INPUT_PATH/*/*_protein.pdb
# rm $INPUT_PATH/*/*_obabel.pdb

# inference:
python inference.py --config=configs_clean/inference.yml
