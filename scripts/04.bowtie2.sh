#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

set -euo pipefail

# ============================================================
# HELIOS NAD-Seq pipeline
# Step 04: Sequential Bowtie2 alignment
#
# Alignment order:
#   1. Spike RNAs
#   2. tRNA + rRNA depletion
#   3. E. coli genome
#
# Usage:
#   sbatch scripts/04.bowtie2.sh \
#       <TRIMMOMATIC_DIR> \
#       <SPIKE_INDEX> \
#       <TRNA_RRNA_INDEX> \
#       <ECOLI_INDEX> \
#       [OUTPUT_DIR]
#
# Example:
#   sbatch scripts/04.bowtie2.sh \
#       results/trimmomatic \
#       /path/to/spikeRnas \
#       /path/to/trna_rrna \
#       /path/to/NC_00913.3_pUC19c
#
# Optional output directory:
#   results/alignment
#
# If OUTPUT_DIR is omitted, "alignment" is created
# next to the Trimmomatic input directory.
# ============================================================


# ------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------

if [ $# -lt 4 ]; then
    echo "Usage:"
    echo "$0 <TRIMMOMATIC_DIR> <SPIKE_INDEX> <TRNA_RRNA_INDEX> <ECOLI_INDEX> [OUTPUT_DIR]"
    exit 1
fi


TRIM_DIR="${1%/}"
SPIKE_INDEX="$2"
TRNA_RRNA_INDEX="$3"
ECOLI_INDEX="$4"


# ------------------------------------------------------------
# Set output directory
# ------------------------------------------------------------

if [ $# -ge 5 ]; then
    OUTPUT_DIR="$5"
else
    parent_dir="$(dirname "$TRIM_DIR")"
    OUTPUT_DIR="${parent_dir}/alignment"
fi

mkdir -p "$OUTPUT_DIR"


# ------------------------------------------------------------
# Check input directory
# ------------------------------------------------------------

if [ ! -d "$TRIM_DIR" ]; then
    echo "ERROR: Trimmomatic input directory does not exist:"
    echo "$TRIM_DIR"
    exit 1
fi


# ------------------------------------------------------------
# Load Bowtie2
# ------------------------------------------------------------

module load bio/bowtie2/2.4.5


echo "========================================"
echo "HELIOS NAD-Seq: Step 04 - Bowtie2"
echo "========================================"
echo "Trimmomatic input: $TRIM_DIR"
echo "Output directory:  $OUTPUT_DIR"
echo
echo "Spike index:       $SPIKE_INDEX"
echo "tRNA/rRNA index:   $TRNA_RRNA_INDEX"
echo "E. coli index:     $ECOLI_INDEX"
echo "========================================"


# ------------------------------------------------------------
# Find paired R1 files from Step 03
# ------------------------------------------------------------

shopt -s nullglob
R1_FILES=("$TRIM_DIR"/*R1_trimmed_paired.fastq)

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "ERROR: No paired R1 FASTQ files found."
    echo "Expected pattern:"
    echo "*R1_trimmed_paired.fastq"
    exit 1
fi


echo "Samples found: ${#R1_FILES[@]}"


# ============================================================
# Process each sample
# ============================================================

for R1_PATH in "${R1_FILES[@]}"; do

    r1="$(basename "$R1_PATH")"

    echo
    echo "========================================"
    echo "Sample: $r1"
    echo "========================================"


    # --------------------------------------------------------
    # Identify paired R2
    # --------------------------------------------------------

    R2_PATH="${R1_PATH/_R1_trimmed_paired.fastq/_R2_trimmed_paired.fastq}"

    if [ ! -f "$R2_PATH" ]; then
        echo "ERROR: Missing paired R2 file for:"
        echo "$R1_PATH"
        continue
    fi


    # --------------------------------------------------------
    # Identify singleton files from Trimmomatic
    # --------------------------------------------------------

    R1_UNPAIRED="${R1_PATH/_paired.fastq/_unpaired.fastq}"
    R2_UNPAIRED="${R2_PATH/_paired.fastq/_unpaired.fastq}"


    if [ ! -f "$R1_UNPAIRED" ]; then
        echo "WARNING: Missing R1 unpaired file:"
        echo "$R1_UNPAIRED"
        : > "$R1_UNPAIRED"
    fi

    if [ ! -f "$R2_UNPAIRED" ]; then
        echo "WARNING: Missing R2 unpaired file:"
        echo "$R2_UNPAIRED"
        : > "$R2_UNPAIRED"
    fi


    # --------------------------------------------------------
    # Sample base name
    # --------------------------------------------------------

    base="${r1%_R1_trimmed_paired.fastq}"


    # Create a separate directory for each sample
    sample_dir="${OUTPUT_DIR}/${base}"
    mkdir -p "$sample_dir"


    # ========================================================
    # STEP 1
    # Align to spike RNA
    # ========================================================

    spike_paired_sam="${sample_dir}/${base}_spike_paired.sam"

    spike_unpaired_r1_sam="${sample_dir}/${base}_spike_unpaired_R1.sam"
    spike_unpaired_r2_sam="${sample_dir}/${base}_spike_unpaired_R2.sam"

    spike_unconc_pref="${sample_dir}/${base}_unaligned_spike_paired"

    spike_R1_unp_fastq="${sample_dir}/${base}_unaligned_spike_R1_unpaired.fastq"
    spike_R2_unp_fastq="${sample_dir}/${base}_unaligned_spike_R2_unpaired.fastq"


    if [[ -f "$spike_paired_sam" &&
          -f "$spike_unpaired_r1_sam" &&
          -f "$spike_unpaired_r2_sam" ]]; then

        echo "Spike alignment already exists. Skipping."

    else

        echo "-> Aligning paired reads to spike RNA"

        bowtie2 \
            -x "$SPIKE_INDEX" \
            -1 "$R1_PATH" \
            -2 "$R2_PATH" \
            -S "$spike_paired_sam" \
            --un-conc "${spike_unconc_pref}.fastq"


        echo "-> Aligning R1 singletons to spike RNA"

        if [ -s "$R1_UNPAIRED" ]; then

            bowtie2 \
                -x "$SPIKE_INDEX" \
                -U "$R1_UNPAIRED" \
                -S "$spike_unpaired_r1_sam" \
                --un "$spike_R1_unp_fastq"

        else

            : > "$spike_unpaired_r1_sam"
            : > "$spike_R1_unp_fastq"

        fi


        echo "-> Aligning R2 singletons to spike RNA"

        if [ -s "$R2_UNPAIRED" ]; then

            bowtie2 \
                -x "$SPIKE_INDEX" \
                -U "$R2_UNPAIRED" \
                -S "$spike_unpaired_r2_sam" \
                --un "$spike_R2_unp_fastq"

        else

            : > "$spike_unpaired_r2_sam"
            : > "$spike_R2_unp_fastq"

        fi

    fi


    # ========================================================
    # STEP 2
    # tRNA + rRNA depletion
    # ========================================================

    rrna_paired_sam="${sample_dir}/${base}_rrna_paired.sam"

    rrna_unpaired_r1_sam="${sample_dir}/${base}_rrna_unpaired_R1.sam"
    rrna_unpaired_r2_sam="${sample_dir}/${base}_rrna_unpaired_R2.sam"

    rrna_unconc_pref="${sample_dir}/${base}_unaligned_rrna_paired"

    rrna_R1_unp_fastq="${sample_dir}/${base}_unaligned_rrna_R1_unpaired.fastq"
    rrna_R2_unp_fastq="${sample_dir}/${base}_unaligned_rrna_R2_unpaired.fastq"

    spike_paired1="${spike_unconc_pref}.1.fastq"
    spike_paired2="${spike_unconc_pref}.2.fastq"


    if [[ -f "$rrna_paired_sam" &&
          -f "$rrna_unpaired_r1_sam" &&
          -f "$rrna_unpaired_r2_sam" ]]; then

        echo "tRNA/rRNA depletion already exists. Skipping."

    else

        echo "-> Aligning spike-unaligned paired reads to tRNA/rRNA"

        if [[ -s "$spike_paired1" &&
              -s "$spike_paired2" ]]; then

            bowtie2 \
                -x "$TRNA_RRNA_INDEX" \
                -1 "$spike_paired1" \
                -2 "$spike_paired2" \
                -S "$rrna_paired_sam" \
                --un-conc "${rrna_unconc_pref}.fastq"

        else

            : > "$rrna_paired_sam"
            : > "${rrna_unconc_pref}.1.fastq"
            : > "${rrna_unconc_pref}.2.fastq"

        fi


        echo "-> Aligning R1 singletons to tRNA/rRNA"

        if [ -s "$spike_R1_unp_fastq" ]; then

            bowtie2 \
                -x "$TRNA_RRNA_INDEX" \
                -U "$spike_R1_unp_fastq" \
                -S "$rrna_unpaired_r1_sam" \
                --un "$rrna_R1_unp_fastq"

        else

            : > "$rrna_unpaired_r1_sam"
            : > "$rrna_R1_unp_fastq"

        fi


        echo "-> Aligning R2 singletons to tRNA/rRNA"

        if [ -s "$spike_R2_unp_fastq" ]; then

            bowtie2 \
                -x "$TRNA_RRNA_INDEX" \
                -U "$spike_R2_unp_fastq" \
                -S "$rrna_unpaired_r2_sam" \
                --un "$rrna_R2_unp_fastq"

        else

            : > "$rrna_unpaired_r2_sam"
            : > "$rrna_R2_unp_fastq"

        fi

    fi


    # ========================================================
    # STEP 3
    # Align remaining reads to E. coli
    # ========================================================

    ecoli_paired_sam="${sample_dir}/${base}_eColi_paired.sam"

    ecoli_R1_unp_sam="${sample_dir}/${base}_eColi_R1_unpaired.sam"
    ecoli_R2_unp_sam="${sample_dir}/${base}_eColi_R2_unpaired.sam"

    ecoli_paired1="${rrna_unconc_pref}.1.fastq"
    ecoli_paired2="${rrna_unconc_pref}.2.fastq"


    if [[ -f "$ecoli_paired_sam" &&
          -f "$ecoli_R1_unp_sam" &&
          -f "$ecoli_R2_unp_sam" ]]; then

        echo "E. coli alignment already exists. Skipping."

    else

        echo "-> Aligning paired reads to E. coli"

        if [[ -s "$ecoli_paired1" &&
              -s "$ecoli_paired2" ]]; then

            bowtie2 \
                -x "$ECOLI_INDEX" \
                -1 "$ecoli_paired1" \
                -2 "$ecoli_paired2" \
                -S "$ecoli_paired_sam" \
                --local

        else

            echo "No paired reads remain after tRNA/rRNA depletion."
            : > "$ecoli_paired_sam"

        fi


        echo "-> Aligning R1 singletons to E. coli"

        if [ -s "$rrna_R1_unp_fastq" ]; then

            bowtie2 \
                -x "$ECOLI_INDEX" \
                -U "$rrna_R1_unp_fastq" \
                -S "$ecoli_R1_unp_sam" \
                --local

        else

            : > "$ecoli_R1_unp_sam"

        fi


        echo "-> Aligning R2 singletons to E. coli"

        if [ -s "$rrna_R2_unp_fastq" ]; then

            bowtie2 \
                -x "$ECOLI_INDEX" \
                -U "$rrna_R2_unp_fastq" \
                -S "$ecoli_R2_unp_sam" \
                --local

        else

            : > "$ecoli_R2_unp_sam"

        fi

    fi


    echo "=== Finished $base ==="

done


echo
echo "========================================"
echo "Step 04 completed."
echo "Output directory: $OUTPUT_DIR"
echo "========================================"
