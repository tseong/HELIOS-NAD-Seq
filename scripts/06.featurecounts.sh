#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

set -euo pipefail

# ============================================================
# HELIOS NAD-Seq pipeline
# Step 06: featureCounts
#
# Counting strategy:
#
#   1. Intergenic annotation
#      Input: A-start filtered E. coli SAM files
#
#   2. Variable 3' annotation
#      Input: original E. coli alignment SAM files from Step 04
#
#
# Usage:
#
#   sbatch scripts/06.featurecounts.sh \
#       <ASTART_DIR> \
#       <ALIGNMENT_DIR> \
#       <INTERGENIC_GTF> \
#       <VARIABLE_3PRIME_GTF> \
#       [OUTPUT_DIR]
#
#
# Example:
#
#   sbatch scripts/06.featurecounts.sh \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/Astart \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/alignment \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/gtf/GCF_000005845.2_genomic_intergenic.gtf \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/gtf/GCF_000005845.2_genomic_variable_3prime.gtf
#
#
# If OUTPUT_DIR is omitted:
#
#   <parent_of_ASTART_DIR>/featurecounts
#
#
# Output:
#
# featurecounts/
# ├── intergenic/
# │   └── *.Astart.table
# │
# └── variable_3prime/
#     └── *.table
#
# ============================================================


# ------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------

if [ $# -lt 4 ]; then
    echo "Usage:"
    echo "$0 <ASTART_DIR> <ALIGNMENT_DIR> <INTERGENIC_GTF> <VARIABLE_3PRIME_GTF> [OUTPUT_DIR]"
    exit 1
fi


ASTART_DIR="${1%/}"
ALIGNMENT_DIR="${2%/}"
INTERGENIC_GTF="$3"
VARIABLE_3PRIME_GTF="$4"


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

if [ $# -ge 5 ]; then
    OUTPUT_DIR="${5%/}"
else
    parent_dir="$(dirname "$ASTART_DIR")"
    OUTPUT_DIR="${parent_dir}/featurecounts"
fi


INTERGENIC_OUT="${OUTPUT_DIR}/intergenic"
VARIABLE3_OUT="${OUTPUT_DIR}/variable_3prime"


mkdir -p "$INTERGENIC_OUT"
mkdir -p "$VARIABLE3_OUT"


# ------------------------------------------------------------
# Check input directories/files
# ------------------------------------------------------------

if [ ! -d "$ASTART_DIR" ]; then
    echo "ERROR: A-start directory does not exist:"
    echo "$ASTART_DIR"
    exit 1
fi


if [ ! -d "$ALIGNMENT_DIR" ]; then
    echo "ERROR: Alignment directory does not exist:"
    echo "$ALIGNMENT_DIR"
    exit 1
fi


if [ ! -f "$INTERGENIC_GTF" ]; then
    echo "ERROR: Intergenic GTF does not exist:"
    echo "$INTERGENIC_GTF"
    exit 1
fi


if [ ! -f "$VARIABLE_3PRIME_GTF" ]; then
    echo "ERROR: Variable 3' GTF does not exist:"
    echo "$VARIABLE_3PRIME_GTF"
    exit 1
fi


# ------------------------------------------------------------
# Locate featureCounts
# ------------------------------------------------------------

FEATURECOUNTS="$(command -v featureCounts || true)"


if [ -z "$FEATURECOUNTS" ]; then
    echo "ERROR: featureCounts is not available."
    echo
    echo "Activate the HELIOS environment:"
    echo
    echo "conda activate /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/eColiHelios_2"
    exit 1
fi


echo "========================================"
echo "HELIOS NAD-Seq: Step 06 - featureCounts"
echo "========================================"
echo
echo "A-start input:"
echo "$ASTART_DIR"
echo
echo "Alignment input:"
echo "$ALIGNMENT_DIR"
echo
echo "Intergenic GTF:"
echo "$INTERGENIC_GTF"
echo
echo "Variable 3' GTF:"
echo "$VARIABLE_3PRIME_GTF"
echo
echo "featureCounts:"
echo "$FEATURECOUNTS"
echo
echo "Output directory:"
echo "$OUTPUT_DIR"
echo
echo "========================================"


# ============================================================
# Helper function
# ============================================================

run_featurecounts() {

    SAM_FILE="$1"
    GTF_FILE="$2"
    FEATURE_TYPE="$3"
    OUTPUT_FILE="$4"
    PAIRED="$5"


    if [ -s "$OUTPUT_FILE" ]; then
        echo "Output already exists:"
        echo "$OUTPUT_FILE"
        echo "Skipping."
        return
    fi


    rm -f "$OUTPUT_FILE"
    rm -f "${OUTPUT_FILE}.summary"


    if [ "$PAIRED" = "yes" ]; then

        "$FEATURECOUNTS" \
            -p \
            -a "$GTF_FILE" \
            -t "$FEATURE_TYPE" \
            -g gene_id \
            -o "$OUTPUT_FILE" \
            "$SAM_FILE"

    else

        "$FEATURECOUNTS" \
            -a "$GTF_FILE" \
            -t "$FEATURE_TYPE" \
            -g gene_id \
            -o "$OUTPUT_FILE" \
            "$SAM_FILE"

    fi
}


# ============================================================
# PART 1
# Intergenic annotation using A-start SAM files
# ============================================================

echo
echo "========================================"
echo "PART 1: Intergenic counts from A-start"
echo "========================================"


shopt -s nullglob

ASTART_SAM_FILES=(
    "$ASTART_DIR"/*_eColi_*.Astart.sam
)


if [ ${#ASTART_SAM_FILES[@]} -eq 0 ]; then

    echo "ERROR: No A-start E. coli SAM files found."
    echo "Expected:"
    echo "*_eColi_*.Astart.sam"
    exit 1

fi


echo "A-start SAM files found: ${#ASTART_SAM_FILES[@]}"


for SAM_FILE in "${ASTART_SAM_FILES[@]}"; do

    filename="$(basename "$SAM_FILE")"
    base="${filename%.sam}"

    echo
    echo "----------------------------------------"
    echo "Intergenic counting:"
    echo "$filename"


    if [[ "$filename" == *_eColi_paired.Astart.sam ]]; then

        paired="yes"

    else

        paired="no"

    fi


    OUTPUT_FILE="${INTERGENIC_OUT}/${base}.table"


    run_featurecounts \
        "$SAM_FILE" \
        "$INTERGENIC_GTF" \
        "intergenic" \
        "$OUTPUT_FILE" \
        "$paired"

done


# ============================================================
# PART 2
# Variable 3' annotation using original alignment SAM files
# ============================================================

echo
echo "========================================"
echo "PART 2: Variable 3' counts from alignment"
echo "========================================"


# ------------------------------------------------------------
# Step 04 alignment files are stored inside sample
# subdirectories, so use find recursively.
# ------------------------------------------------------------

mapfile -t ALIGNMENT_SAM_FILES < <(
    find "$ALIGNMENT_DIR" \
        -type f \
        \( \
            -name "*_eColi_paired.sam" \
            -o \
            -name "*_eColi_R1_unpaired.sam" \
            -o \
            -name "*_eColi_R2_unpaired.sam" \
        \) \
        | sort
)


if [ ${#ALIGNMENT_SAM_FILES[@]} -eq 0 ]; then

    echo "ERROR: No E. coli SAM files found recursively under:"
    echo "$ALIGNMENT_DIR"
    exit 1

fi


echo "Alignment SAM files found: ${#ALIGNMENT_SAM_FILES[@]}"


for SAM_FILE in "${ALIGNMENT_SAM_FILES[@]}"; do

    filename="$(basename "$SAM_FILE")"
    base="${filename%.sam}"


    echo
    echo "----------------------------------------"
    echo "Variable 3' counting:"
    echo "$filename"


    if [[ "$filename" == *_eColi_paired.sam ]]; then

        paired="yes"

    else

        paired="no"

    fi


    OUTPUT_FILE="${VARIABLE3_OUT}/${base}.table"


    run_featurecounts \
        "$SAM_FILE" \
        "$VARIABLE_3PRIME_GTF" \
        "gene" \
        "$OUTPUT_FILE" \
        "$paired"

done


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "========================================"
echo "Step 06 completed."
echo "========================================"
echo
echo "Intergenic counts:"
echo "$INTERGENIC_OUT"
echo
echo "Variable 3' counts:"
echo "$VARIABLE3_OUT"
echo
echo "========================================"
