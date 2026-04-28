#!/bin/bash

# Prepares targets for DiffDock.

# Delete stuff to start fresh
rm /home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_lig/*.sdf
rm /home/daniil/Desktop/Ph.D/DiffDock/data/targets_*.csv
rm -rf /home/daniil/Desktop/Ph.D/DiffDock/results/user_inference/*

# data preprocessing:
python prepare_targets.py
echo "Number of total converted molecules:"
ls /home/daniil/Desktop/Ph.D/ThesisRepo/databases/QBioLiP/nr_lig/*.sdf | wc -l

# inference:
cd /home/daniil/Desktop/Ph.D/DiffDock/

# the prepare_targets.py file creates small csv files which allows to load data in batches
CSVS=/home/daniil/Desktop/Ph.D/DiffDock/data/targets_*.csv
for CSV in $CSVS
do
	python -m inference --config default_inference_args.yaml --protein_ligand_csv $CSV
done