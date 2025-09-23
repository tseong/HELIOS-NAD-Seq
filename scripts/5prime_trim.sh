#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

#trim the 5' end of the reads
source /opt/bwhpc/common/devel/miniconda/3-py39-4.12.0/etc/profile.d/conda.sh
conda activate env.helios.yml
python scripts/5prime_trim.py
