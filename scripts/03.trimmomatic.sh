#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

set -euo pipefail

# ============================================================
# HELIOS NAD-Seq pipeline
# Step 03: Trimmomatic paired-end trimming
#
# Usage:
#   sbatch scripts/03.trimmomatic.sh <INPUT_DIR> [OUTPUT_DIR]
#
# Example:
#   sbatch scripts/03.trimmomatic.sh \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/5prime_trimmed
#
# If OUTPUT_DIR is omitted, "trimmomatic" will be created
# next to the input directory.
#
# Multiple copies of this script may run concurrently.
# A per-sample lock directory prevents duplicate processing.
# ============================================================


# ------------------------------------------------------------
# Check input arguments
# ------------------------------------------------------------

if [ $# -lt 1 ]; then
    echo "Usage: $0 <INPUT_DIR> [OUTPUT_DIR]"
    exit 1
fi


# ------------------------------------------------------------
# Input directory
# ------------------------------------------------------------

input_dir="${1%/}"

if [ ! -d "$input_dir" ]; then
    echo "ERROR: Input directory does not exist:"
    echo "$input_dir"
    exit 1
fi


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

if [ $# -ge 2 ]; then
    output_dir="${2%/}"
else
    parent_dir="$(dirname "$input_dir")"
    output_dir="${parent_dir}/trimmomatic"
fi

mkdir -p "$output_dir"


# ------------------------------------------------------------
# Check active Conda environment
# ------------------------------------------------------------

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: No Conda environment is active."
    exit 1
fi


# ------------------------------------------------------------
# Locate Trimmomatic
# ------------------------------------------------------------

TRIMMOMATIC="$(command -v trimmomatic || true)"

if [ -z "$TRIMMOMATIC" ]; then
    echo "ERROR: Trimmomatic is not available."
    echo "Active Conda environment:"
    echo "$CONDA_PREFIX"
    exit 1
fi


# ------------------------------------------------------------
# Locate TruSeq3-PE adapter file
# ------------------------------------------------------------

ADAPTERS="$(find "$CONDA_PREFIX" \
    -name "TruSeq3-PE.fa" \
    -type f \
    2>/dev/null \
    | head -n 1)"

if [ -z "$ADAPTERS" ] || [ ! -f "$ADAPTERS" ]; then
    echo "ERROR: TruSeq3-PE.fa could not be found."
    echo "Searched in:"
    echo "$CONDA_PREFIX"
    exit 1
fi


faPara="1:30:15"


# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

log_dir="${output_dir}/logs"
lock_dir="${output_dir}/locks"

mkdir -p "$log_dir"
mkdir -p "$lock_dir"


echo "========================================"
echo "HELIOS NAD-Seq: Step 03 - Trimmomatic"
echo "========================================"
echo "Input directory:   $input_dir"
echo "Output directory:  $output_dir"
echo "Lock directory:    $lock_dir"
echo "Conda environment: $CONDA_PREFIX"
echo "Trimmomatic:       $TRIMMOMATIC"
echo "Adapters:          $ADAPTERS"
echo "========================================"


# ------------------------------------------------------------
# Find Step 02 R1 files
# ------------------------------------------------------------

shopt -s nullglob

R1_files=(
    "$input_dir"/*R1_trimmed.fastq
)

if [ ${#R1_files[@]} -eq 0 ]; then
    echo "ERROR: No R1 trimmed FASTQ files found in:"
    echo "$input_dir"
    echo
    echo "Expected pattern:"
    echo "*R1_trimmed.fastq"
    exit 1
fi


echo "R1 FASTQ files found: ${#R1_files[@]}"
echo


# ============================================================
# Process each sample
# ============================================================

for R1_file in "${R1_files[@]}"; do

    file="$(basename "$R1_file")"

    R2_file="${R1_file/R1_trimmed.fastq/R2_trimmed.fastq}"

    if [ ! -f "$R2_file" ]; then
        echo "WARNING: Matching R2 file not found:"
        echo "$R1_file"
        echo "Expected:"
        echo "$R2_file"
        echo "Skipping..."
        continue
    fi


    # --------------------------------------------------------
    # Sample name
    # --------------------------------------------------------

    base_name="${file%_R1_trimmed.fastq}"


    # --------------------------------------------------------
    # Output files
    # --------------------------------------------------------

    R1_paired="${output_dir}/${base_name}_R1_trimmed_paired.fastq"
    R1_unpaired="${output_dir}/${base_name}_R1_trimmed_unpaired.fastq"

    R2_paired="${output_dir}/${base_name}_R2_trimmed_paired.fastq"
    R2_unpaired="${output_dir}/${base_name}_R2_trimmed_unpaired.fastq"

    trimlog="${log_dir}/${base_name}_trimmomatic.trimlog"


    # --------------------------------------------------------
    # Skip samples already completed
    # --------------------------------------------------------

    if [[ -s "$R1_paired" &&
          -s "$R1_unpaired" &&
          -s "$R2_paired" &&
          -s "$R2_unpaired" ]]; then

        echo "----------------------------------------"
        echo "Already completed: $base_name"
        echo "Skipping..."

        continue
    fi


    # --------------------------------------------------------
    # Attempt to claim sample using atomic mkdir
    #
    # mkdir succeeds for only one concurrent process.
    # All other jobs will see the existing lock and skip.
    # --------------------------------------------------------

    sample_lock="${lock_dir}/${base_name}.lock"

    if ! mkdir "$sample_lock" 2>/dev/null; then

        echo "----------------------------------------"
        echo "Another job is currently processing:"
        echo "$base_name"
        echo "Skipping..."

        continue
    fi


    # --------------------------------------------------------
    # We now own this sample.
    #
    # Remove the lock automatically if this job exits
    # unexpectedly while processing this sample.
    # --------------------------------------------------------

    cleanup_lock() {
        rm -rf "$sample_lock"
    }

    trap cleanup_lock EXIT INT TERM


    # --------------------------------------------------------
    # Re-check outputs after acquiring lock
    #
    # Another job could have completed the sample between
    # the first check and acquisition of the lock.
    # --------------------------------------------------------

    if [[ -s "$R1_paired" &&
          -s "$R1_unpaired" &&
          -s "$R2_paired" &&
          -s "$R2_unpaired" ]]; then

        echo "----------------------------------------"
        echo "Sample completed by another job:"
        echo "$base_name"
        echo "Skipping..."

        rm -rf "$sample_lock"
        trap - EXIT INT TERM

        continue
    fi


    echo "----------------------------------------"
    echo "Processing: $base_name"
    echo "R1: $R1_file"
    echo "R2: $R2_file"
    echo "Lock acquired: $sample_lock"


    # --------------------------------------------------------
    # Remove incomplete outputs from an earlier failed run
    # --------------------------------------------------------

    rm -f \
        "$R1_paired" \
        "$R1_unpaired" \
        "$R2_paired" \
        "$R2_unpaired" \
        "$trimlog"


    # --------------------------------------------------------
    # Run Trimmomatic
    # --------------------------------------------------------

    "$TRIMMOMATIC" \
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


    # --------------------------------------------------------
    # Verify expected output files
    # --------------------------------------------------------

    if [[ ! -f "$R1_paired" ||
          ! -f "$R1_unpaired" ||
          ! -f "$R2_paired" ||
          ! -f "$R2_unpaired" ]]; then

        echo "ERROR: Trimmomatic did not generate all expected files for:"
        echo "$base_name"

        exit 1
    fi


    echo "Completed: $base_name"


    # --------------------------------------------------------
    # Release lock
    # --------------------------------------------------------

    rm -rf "$sample_lock"
    trap - EXIT INT TERM

    echo "Lock released."
    echo

done


echo "========================================"
echo "Step 03 job completed."
echo "Output directory:"
echo "$output_dir"
echo "========================================"
