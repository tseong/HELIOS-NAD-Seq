#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

# ============================================================
# HELIOS NAD-Seq pipeline
# Step 02: 5' trimming
#
# Usage:
#   sbatch scripts/02_5prime_trim_slurm.sh <CUTADAPT_DIR> [OUTPUT_DIR]
#
# Example:
#   sbatch scripts/02_5prime_trim_slurm.sh results/cutadapt
#
# Optional output directory:
#   sbatch scripts/02_5prime_trim_slurm.sh \
#       results/cutadapt \
#       results/5prime_trimmed
# ============================================================


# Check that an input directory was provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <CUTADAPT_DIR> [OUTPUT_DIR]"
    exit 1
fi


# Input and optional output directories
input_dir="$1"
output_dir="$2"


# Activate Conda environment
source /opt/bwhpc/common/devel/miniconda/3-py39-4.12.0/etc/profile.d/conda.sh
conda activate env.helios.yml


echo "========================================"
echo "HELIOS NAD-Seq: Step 02 - 5' trimming"
echo "========================================"
echo "Input directory:  $input_dir"

if [ -n "$output_dir" ]; then
    echo "Output directory: $output_dir"
else
    echo "Output directory: automatically generated"
fi

echo "========================================"


# Run 5' trimming
if [ -n "$output_dir" ]; then

    python scripts/02_5prime_trim.py \
        "$input_dir" \
        "$output_dir"

else

    python scripts/02_5prime_trim.py \
        "$input_dir"

fi


echo "========================================"
echo "Step 02 completed."
echo "========================================"
