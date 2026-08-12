#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

# ============================================================
# Usage:
#   sbatch 02_5prime_trim.sh <CUTADAPT_DIR> [OUTPUT_DIR]
#
# Example:
#   sbatch 02_5prime_trim.sh /path/to/cutadapt
#
# or:
#   sbatch 02_5prime_trim.sh /path/to/cutadapt /path/to/5prime_trimmed
# ============================================================

if [ $# -lt 1 ]; then
    echo "Usage: $0 <CUTADAPT_DIR> [OUTPUT_DIR]"
    exit 1
fi

input_dir="$1"
output_dir="$2"

source /opt/bwhpc/common/devel/miniconda/3-py39-4.12.0/etc/profile.d/conda.sh
conda activate env.helios.yml

if [ -n "$output_dir" ]; then
    python scripts/5prime_trim.py "$input_dir" "$output_dir"
else
    python scripts/5prime_trim.py "$input_dir"
fi
