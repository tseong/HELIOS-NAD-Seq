#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

set -euo pipefail

# ============================================================
# Download all SRA runs from study SRP619047
#
# Output:
#   /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/fastq/
#
# Requires:
#   prefetch
#   fasterq-dump
#   esearch
#   efetch
#
# The script:
#   1. Retrieves all SRR accessions associated with SRP619047
#   2. Skips runs whose FASTQ files already exist
#   3. Downloads missing runs with prefetch
#   4. Converts each run to FASTQ
#   5. Renames paired files to:
#        SRRxxxx_R1_001.fastq
#        SRRxxxx_R2_001.fastq
# ============================================================


WORKDIR="/gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow"

SRA_DIR="${WORKDIR}/sra"
FASTQ_DIR="${WORKDIR}/fastq"
TMP_DIR="${WORKDIR}/tmp"

RUN_LIST="${WORKDIR}/SRP619047_runs.txt"

STUDY="SRP619047"


# ------------------------------------------------------------
# Create directories
# ------------------------------------------------------------

mkdir -p "$SRA_DIR"
mkdir -p "$FASTQ_DIR"
mkdir -p "$TMP_DIR"


# ------------------------------------------------------------
# Check required programs
# ------------------------------------------------------------

if ! command -v prefetch >/dev/null 2>&1; then
    echo "ERROR: prefetch is not available."
    echo "Activate an environment containing sra-tools."
    exit 1
fi

if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "ERROR: fasterq-dump is not available."
    echo "Activate an environment containing sra-tools."
    exit 1
fi

if ! command -v esearch >/dev/null 2>&1; then
    echo "ERROR: esearch is not available."
    echo "Install NCBI Entrez Direct (entrez-direct)."
    exit 1
fi

if ! command -v efetch >/dev/null 2>&1; then
    echo "ERROR: efetch is not available."
    echo "Install NCBI Entrez Direct (entrez-direct)."
    exit 1
fi


echo "========================================"
echo "NCBI SRA project download"
echo "========================================"
echo "Study:      $STUDY"
echo "Workspace:  $WORKDIR"
echo "========================================"


# ------------------------------------------------------------
# Step 1: Retrieve all SRR accessions
# ------------------------------------------------------------

echo
echo "Retrieving run accessions for $STUDY..."

esearch \
    -db sra \
    -query "$STUDY" \
    | efetch \
        -format runinfo \
    | cut -d',' -f1 \
    | grep '^SRR' \
    | sort -u \
    > "$RUN_LIST"


RUN_COUNT=$(wc -l < "$RUN_LIST")


if [ "$RUN_COUNT" -eq 0 ]; then
    echo "ERROR: No SRR accessions found for $STUDY."
    exit 1
fi


echo "Found $RUN_COUNT runs:"
cat "$RUN_LIST"


# ============================================================
# Step 2: Download and convert each run
# ============================================================

while read -r RUN; do

    echo
    echo "========================================"
    echo "Processing $RUN"
    echo "========================================"


    NEW_R1="${FASTQ_DIR}/${RUN}_R1_001.fastq"
    NEW_R2="${FASTQ_DIR}/${RUN}_R2_001.fastq"
    NEW_SINGLE="${FASTQ_DIR}/${RUN}_R1_001.fastq"

    RAW_R1="${FASTQ_DIR}/${RUN}_1.fastq"
    RAW_R2="${FASTQ_DIR}/${RUN}_2.fastq"
    RAW_SINGLE="${FASTQ_DIR}/${RUN}.fastq"


    # --------------------------------------------------------
    # Skip already completed runs
    # --------------------------------------------------------

    if [[ -s "$NEW_R1" && -s "$NEW_R2" ]]; then

        echo "Paired FASTQ files already exist."
        echo "Skipping $RUN."

        continue

    elif [[ -s "$NEW_SINGLE" && ! -e "$NEW_R2" ]]; then

        echo "Single-end FASTQ file already exists."
        echo "Skipping $RUN."

        continue
    fi


    # --------------------------------------------------------
    # Handle partially converted FASTQ output
    # --------------------------------------------------------

    if [[ -s "$RAW_R1" && -s "$RAW_R2" ]]; then

        echo "Unrenamed paired FASTQ files already exist."
        echo "Renaming without downloading again."

        mv "$RAW_R1" "$NEW_R1"
        mv "$RAW_R2" "$NEW_R2"

        continue

    elif [[ -s "$RAW_SINGLE" ]]; then

        echo "Unrenamed single-end FASTQ already exists."
        echo "Renaming without downloading again."

        mv "$RAW_SINGLE" "$NEW_SINGLE"

        continue
    fi


    # --------------------------------------------------------
    # Download SRA archive only if it is not already present
    # --------------------------------------------------------

    SRA_RUN_DIR="${SRA_DIR}/${RUN}"

    if [ -d "$SRA_RUN_DIR" ]; then

        echo "SRA archive already exists for $RUN."
        echo "Skipping prefetch."

    else

        echo "Downloading $RUN..."

        prefetch "$RUN" \
            --output-directory "$SRA_DIR"

    fi


    # --------------------------------------------------------
    # Convert SRA to FASTQ
    # --------------------------------------------------------

    echo "Converting $RUN to FASTQ..."

    fasterq-dump \
        "$SRA_RUN_DIR" \
        --split-files \
        --threads 8 \
        --temp "$TMP_DIR" \
        --outdir "$FASTQ_DIR"


    # --------------------------------------------------------
    # Rename paired-end or single-end output
    # --------------------------------------------------------

    if [[ -s "$RAW_R1" && -s "$RAW_R2" ]]; then

        mv "$RAW_R1" "$NEW_R1"
        mv "$RAW_R2" "$NEW_R2"

        echo "Paired-end FASTQ created:"
        echo "  $NEW_R1"
        echo "  $NEW_R2"

    elif [ -s "$RAW_SINGLE" ]; then

        mv "$RAW_SINGLE" "$NEW_SINGLE"

        echo "Single-end FASTQ created:"
        echo "  $NEW_SINGLE"

    else

        echo "ERROR: Expected FASTQ output not found for $RUN."
        exit 1
    fi


done < "$RUN_LIST"


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "========================================"
echo "Project download completed."
echo "========================================"
echo "Study:          $STUDY"
echo "Runs listed:    $RUN_COUNT"
echo "FASTQ directory:"
echo "$FASTQ_DIR"
echo
echo "Downloaded FASTQ files:"

ls -lh "$FASTQ_DIR"
