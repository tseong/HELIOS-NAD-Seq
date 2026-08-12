#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

# ============================================================
# HELIOS NAD-Seq pipeline
# Step 03: Trimmomatic paired-end trimming
#
# Usage:
#   sbatch scripts/03_trimmomatic.sh <INPUT_DIR> [OUTPUT_DIR]
#
# Example:
#   sbatch scripts/03_trimmomatic.sh results/5prime_trimmed
#
# Optional output directory:
#   sbatch scripts/03_trimmomatic.sh \
#       results/5prime_trimmed \
#       results/trimmomatic
#
# If OUTPUT_DIR is not provided, "trimmomatic" will be created
# next to the input directory.
# ============================================================


# Check input argument
if [ $# -lt 1 ]; then
    echo "Usage: $0 <INPUT_DIR> [OUTPUT_DIR]"
    exit 1
fi


# Input directory from Step 02
input_dir="$1"
input_dir="${input_dir%/}"


# Check input directory
if [ ! -d "$input_dir" ]; then
    echo "ERROR: Input directory does not exist:"
    echo "$input_dir"
    exit 1
fi


# Set output directory
if [ $# -ge 2 ]; then
    output_dir="$2"
else
    parent_dir="$(dirname "$input_dir")"
    output_dir="${parent_dir}/trimmomatic"
fi


mkdir -p "$output_dir"


echo "========================================"
echo "HELIOS NAD-Seq: Step 03 - Trimmomatic"
echo "========================================"
echo "Input directory:  $input_dir"
echo "Output directory: $output_dir"
echo "========================================"


# ============================================================
# Trimmomatic configuration
#
# Modify these paths according to your installation.
# ============================================================

JAVA="/home/hd/hd_hd/hd_uv268/software/jdk1.8.0_361/bin/java"

TRIMMOMATIC_JAR="/home/hd/hd_hd/hd_uv268/software/Trimmomatic-0.39/trimmomatic-0.39.jar"

ADAPTERS="/home/hd/hd_hd/hd_uv268/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

faPara="1:30:15"


# Check required files
if [ ! -x "$JAVA" ]; then
    echo "ERROR: Java executable not found:"
    echo "$JAVA"
    exit 1
fi

if [ ! -f "$TRIMMOMATIC_JAR" ]; then
    echo "ERROR: Trimmomatic JAR not found:"
    echo "$TRIMMOMATIC_JAR"
    exit 1
fi

if [ ! -f "$ADAPTERS" ]; then
    echo "ERROR: Trimmomatic adapter file not found:"
    echo "$ADAPTERS"
    exit 1
fi


# Log directory
log_dir="${output_dir}/logs"
mkdir -p "$log_dir"


# Find Step 02 R1 files
shopt -s nullglob
R1_files=("$input_dir"/*R1_001_trimmed.fastq)

if [ ${#R1_files[@]} -eq 0 ]; then
    echo "ERROR: No R1 trimmed FASTQ files found in:"
    echo "$input_dir"
    echo "Expected pattern: *R1_001_trimmed.fastq"
    exit 1
fi


echo "FASTQ pairs found: ${#R1_files[@]}"


for R1_file in "${R1_files[@]}"; do

    # File name only
    file="$(basename "$R1_file")"

    # Identify corresponding R2
    R2_file="${R1_file/R1_001_trimmed.fastq/R2_001_trimmed.fastq}"

    if [ ! -f "$R2_file" ]; then
        echo "WARNING: Matching R2 file not found for:"
        echo "$R1_file"
        echo "Skipping..."
        continue
    fi


    # Sample/base name
    base_name="${file%_R1_001_trimmed.fastq}"


    # Output files
    R1_paired="${output_dir}/${base_name}_R1_trimmed_paired.fastq"
    R1_unpaired="${output_dir}/${base_name}_R1_trimmed_unpaired.fastq"

    R2_paired="${output_dir}/${base_name}_R2_trimmed_paired.fastq"
    R2_unpaired="${output_dir}/${base_name}_R2_trimmed_unpaired.fastq"

    trimlog="${log_dir}/${base_name}_trimmomatic.trimlog"


    # Skip completed samples
    if [[ -f "$R1_paired" &&
          -f "$R1_unpaired" &&
          -f "$R2_paired" &&
          -f "$R2_unpaired" ]]; then

        echo "Trimmomatic output for $base_name already exists. Skipping..."
        continue
    fi


    echo "----------------------------------------"
    echo "Processing: $base_name"
    echo "R1: $R1_file"
    echo "R2: $R2_file"


    "$JAVA" \
        -jar "$TRIMMOMATIC_JAR" \
        PE \
        -phred33 \
        -trimlog "$trimlog" \
        "$R1_file" \
        "$R2_file" \
        "$R1_paired" \
        "$R1_unpaired" \
        "$R2_paired" \
        "$R2_unpaired" \
        ILLUMINACLIP:"${ADAPTERS}:${faPara}" \
        SLIDINGWINDOW:4:20 \
        MINLEN:18


    # Check exit status
    if [ $? -ne 0 ]; then
        echo "ERROR: Trimmomatic failed for $base_name"
        exit 1
    fi

done


echo "========================================"
echo "Step 03 completed."
echo "Output directory: $output_dir"
echo "========================================"
