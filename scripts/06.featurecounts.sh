```bash
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
# Usage:
#   sbatch scripts/06.featurecounts.sh \
#       <ASTART_DIR> \
#       <INTERGENIC_GTF> \
#       [OUTPUT_DIR]
#
# Example:
#   sbatch scripts/06.featurecounts.sh \
#       results/Astart \
#       reference/ecoli_intergenic.gtf
#
# If OUTPUT_DIR is omitted, "featurecounts" will be created
# next to the Astart input directory.
# ============================================================


# ------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------

if [ $# -lt 2 ]; then
    echo "Usage:"
    echo "$0 <ASTART_DIR> <INTERGENIC_GTF> [OUTPUT_DIR]"
    exit 1
fi


ASTART_DIR="${1%/}"
INTERGENIC_GTF="$2"


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

if [ $# -ge 3 ]; then
    OUTPUT_DIR="$3"
else
    parent_dir="$(dirname "$ASTART_DIR")"
    OUTPUT_DIR="${parent_dir}/featurecounts"
fi


mkdir -p "$OUTPUT_DIR"


# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

if [ ! -d "$ASTART_DIR" ]; then
    echo "ERROR: A-start SAM directory does not exist:"
    echo "$ASTART_DIR"
    exit 1
fi


if [ ! -f "$INTERGENIC_GTF" ]; then
    echo "ERROR: Intergenic GTF does not exist:"
    echo "$INTERGENIC_GTF"
    exit 1
fi


# ------------------------------------------------------------
# featureCounts executable
#
# Change this path if featureCounts is installed elsewhere.
# ------------------------------------------------------------

FEATURECOUNTS="/home/hd/hd_hd/hd_uv268/software/subread-2.0.6-Linux-x86_64/bin/featureCounts"


if [ ! -x "$FEATURECOUNTS" ]; then
    echo "ERROR: featureCounts executable not found:"
    echo "$FEATURECOUNTS"
    exit 1
fi


echo "========================================"
echo "HELIOS NAD-Seq: Step 06 - featureCounts"
echo "========================================"
echo "Input SAM directory: $ASTART_DIR"
echo "Intergenic GTF:      $INTERGENIC_GTF"
echo "Output directory:    $OUTPUT_DIR"
echo "========================================"


# ------------------------------------------------------------
# Find A-start E. coli SAM files
# ------------------------------------------------------------

shopt -s nullglob

SAM_FILES=("$ASTART_DIR"/*_eColi_*.Astart.sam)


if [ ${#SAM_FILES[@]} -eq 0 ]; then
    echo "ERROR: No A-start E. coli SAM files found."
    echo "Expected pattern:"
    echo "*_eColi_*.Astart.sam"
    exit 1
fi


echo "SAM files found: ${#SAM_FILES[@]}"
echo


# ============================================================
# Run featureCounts
# ============================================================

for SAM_FILE in "${SAM_FILES[@]}"; do

    filename="$(basename "$SAM_FILE")"

    base="${filename%.sam}"

    OUTPUT_FILE="${OUTPUT_DIR}/${base}.table"


    if [ -f "$OUTPUT_FILE" ]; then
        echo "Skipping $filename:"
        echo "$OUTPUT_FILE already exists."
        echo
        continue
    fi


    echo "----------------------------------------"
    echo "Processing: $filename"


    # --------------------------------------------------------
    # Paired-end SAM
    # --------------------------------------------------------

    if [[ "$filename" == *_paired.Astart.sam ]]; then

        echo "Detected paired-end alignment."

        "$FEATURECOUNTS" \
            -p \
            -a "$INTERGENIC_GTF" \
            -t intergenic \
            -o "$OUTPUT_FILE" \
            "$SAM_FILE"


    # --------------------------------------------------------
    # Single-end SAM
    # --------------------------------------------------------

    else

        echo "Detected single-end alignment."

        "$FEATURECOUNTS" \
            -a "$INTERGENIC_GTF" \
            -t intergenic \
            -o "$OUTPUT_FILE" \
            "$SAM_FILE"

    fi


    # --------------------------------------------------------
    # Check whether featureCounts succeeded
    # --------------------------------------------------------

    if [ $? -ne 0 ]; then
        echo "ERROR: featureCounts failed for:"
        echo "$filename"
        exit 1
    fi


    echo "Output: $OUTPUT_FILE"
    echo

done


echo "========================================"
echo "Step 06 completed."
echo "Output directory: $OUTPUT_DIR"
echo "========================================"
```
