# Machine Learning Docking Review Code Repository
This repository contains all necessary scripts and information to install, preprocess a locally downloaded snapshot of Q-BioLiP database.

The following sections explain how to download Q-BioLiP and install different docking tools (Vina and ML-based).

## Prerequisite
### Conda environement and utilities
1. [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) if not already. Then, create a new environment with the right python version and activate it:
    - `conda create --name ml_review python=3.10`
    - `conda activate ml_review`
2. Install last version of torch from [their website](https://pytorch.org/get-started/locally/). GPU acceleration is recommanded but not necessary for quick replication, required for complete replication.
3. Install required packages using pip:
    - `pip install -r requirements.txt`
4. Install usalign:
    - `conda install -c bioconda usalign`
5. Install obabel:
    - `conda install conda-forge::openbabel`

## Quick replication of paper's results
The quick replication doesn't require to install all Vina and ML-based docking tools but only to compute metrics from already generated outputs. For a complete replication with full output generation please refer to "Complete pipeline".

### Database
The database include complexes computed by different docking tools as well as a minimalist version of Q-BioLiP pre-processed dataset.

#### 1. Download
Download the [zip file](https://tubcloud.tu-berlin.de/s/KF9i7qEwEaGx6gS) and extract the content into database folder.

#### 2. Compute metrics
Run `src/results/combine.py`

#### 4. Explore results
Follow jupyter notebook in `src/results/viz.ipynb`