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
#
#   1. Spike RNAs
#   2. E. coli genome
#
# Only reads that remain unaligned to spike RNA are passed
# to the E. coli genome alignment.
#
# Multiple instances of this script may run simultaneously.
# A per-sample lock prevents two jobs from processing the same
# sample at the same time.
#
#
# Usage:
#
# sbatch scripts/04.bowtie2.sh \
#     <TRIMMOMATIC_DIR> \
#     <SPIKE_INDEX> \
#     <ECOLI_INDEX> \
#     [OUTPUT_DIR]
#
#
# Example:
#
# sbatch scripts/04.bowtie2.sh \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/trimmomatic \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/bowtie2_index/spike \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/bowtie2_index/NC_00913.3_pUC19c
#
# ============================================================


# ------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------

if [ $# -lt 3 ]; then

    echo "Usage:"
    echo "$0 <TRIMMOMATIC_DIR> <SPIKE_INDEX> <ECOLI_INDEX> [OUTPUT_DIR]"

    exit 1

fi


TRIM_DIR="${1%/}"
SPIKE_INDEX="$2"
ECOLI_INDEX="$3"


# ------------------------------------------------------------
# Set output directory
# ------------------------------------------------------------

if [ $# -ge 4 ]; then

    OUTPUT_DIR="${4%/}"

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
# Conda environment / Bowtie2
# ------------------------------------------------------------

CONDA_ENV="/gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/eColiHelios_2"

BOWTIE2="$CONDA_ENV/bin/bowtie2"


# ------------------------------------------------------------
# Use system Perl
#
# Conda Perl previously caused libnsl.so.1 problems on Helix.
# ------------------------------------------------------------

if [ -x /usr/bin/perl ]; then

    PERL="/usr/bin/perl"

elif [ -x /bin/perl ]; then

    PERL="/bin/perl"

else

    echo "ERROR: System Perl not found."

    exit 1

fi


# ------------------------------------------------------------
# Runtime library path
# ------------------------------------------------------------

export LD_LIBRARY_PATH="/lib64:/usr/lib64:${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"


# ------------------------------------------------------------
# Check executables
# ------------------------------------------------------------

if [ ! -f "$BOWTIE2" ]; then

    echo "ERROR: Bowtie2 executable not found:"
    echo "$BOWTIE2"

    exit 1

fi


if [ ! -r "$BOWTIE2" ]; then

    echo "ERROR: Bowtie2 executable is not readable:"
    echo "$BOWTIE2"

    exit 1

fi


if [ ! -x "$PERL" ]; then

    echo "ERROR: System Perl not found:"
    echo "$PERL"

    exit 1

fi


# ------------------------------------------------------------
# Check Bowtie2 indices
# ------------------------------------------------------------

check_index() {

    index_prefix="$1"

    if [[ ! -f "${index_prefix}.1.bt2" &&
          ! -f "${index_prefix}.1.bt2l" ]]; then

        echo "ERROR: Bowtie2 index not found for prefix:"
        echo "$index_prefix"

        exit 1

    fi
}


check_index "$SPIKE_INDEX"
check_index "$ECOLI_INDEX"


# ------------------------------------------------------------
# Runtime diagnostics
# ------------------------------------------------------------

echo "========================================"
echo "HELIOS NAD-Seq: Step 04 - Bowtie2"
echo "========================================"

echo "Trimmomatic input:   $TRIM_DIR"
echo "Output directory:    $OUTPUT_DIR"

echo

echo "Conda environment:   $CONDA_ENV"
echo "Bowtie2 executable:  $BOWTIE2"
echo "System Perl:         $PERL"

echo

echo "Spike index:         $SPIKE_INDEX"
echo "E. coli index:       $ECOLI_INDEX"

echo "========================================"

echo
echo "Testing system Perl:"
"$PERL" --version | head -n 2

echo
echo "Testing Bowtie2:"
"$PERL" "$BOWTIE2" --version

echo
echo "========================================"


# ------------------------------------------------------------
# Find paired R1 files
# ------------------------------------------------------------

shopt -s nullglob


R1_FILES=(
    "$TRIM_DIR"/*R1_trimmed_paired.fastq
)


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

    base="${r1%_R1_trimmed_paired.fastq}"


    echo
    echo "========================================"
    echo "Sample: $base"
    echo "========================================"


    # --------------------------------------------------------
    # Identify paired R2
    # --------------------------------------------------------

    R2_PATH="${R1_PATH/_R1_trimmed_paired.fastq/_R2_trimmed_paired.fastq}"


    if [ ! -f "$R2_PATH" ]; then

        echo "ERROR: Missing paired R2 file:"
        echo "$R2_PATH"

        echo "Skipping..."

        continue

    fi


    # --------------------------------------------------------
    # Identify singleton files
    # --------------------------------------------------------

    R1_UNPAIRED="${R1_PATH/_paired.fastq/_unpaired.fastq}"

    R2_UNPAIRED="${R2_PATH/_paired.fastq/_unpaired.fastq}"


    if [ ! -f "$R1_UNPAIRED" ]; then

        echo "WARNING: Missing R1 unpaired file:"
        echo "$R1_UNPAIRED"

        echo "Creating empty placeholder."

        : > "$R1_UNPAIRED"

    fi


    if [ ! -f "$R2_UNPAIRED" ]; then

        echo "WARNING: Missing R2 unpaired file:"
        echo "$R2_UNPAIRED"

        echo "Creating empty placeholder."

        : > "$R2_UNPAIRED"

    fi


    # --------------------------------------------------------
    # Sample-specific output directory
    # --------------------------------------------------------

    sample_dir="${OUTPUT_DIR}/${base}"

    mkdir -p "$sample_dir"


    # ========================================================
    # SAMPLE LOCK
    # ========================================================

    lock_dir="${sample_dir}/.bowtie2.lock"


    if ! mkdir "$lock_dir" 2>/dev/null; then

        echo "LOCKED by another job."
        echo "Skipping sample: $base"

        continue

    fi


    {
        echo "sample=$base"
        echo "job_id=${SLURM_JOB_ID:-unknown}"
        echo "hostname=$(hostname)"
        echo "started=$(date '+%Y-%m-%d %H:%M:%S')"

    } > "${lock_dir}/info.txt"


    echo "Lock acquired: $lock_dir"


    cleanup_lock() {

        if [ -n "${lock_dir:-}" ] && [ -d "$lock_dir" ]; then

            rm -rf "$lock_dir"

        fi

    }


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


    if [[ -s "$spike_paired_sam" &&
          -s "${spike_unconc_pref}.1.fastq" &&
          -s "${spike_unconc_pref}.2.fastq" &&
          -f "$spike_unpaired_r1_sam" &&
          -f "$spike_unpaired_r2_sam" &&
          -f "$spike_R1_unp_fastq" &&
          -f "$spike_R2_unp_fastq" ]]; then

        echo "Spike alignment already complete. Skipping."

    else

        echo "-> Aligning paired reads to spike RNA"


        "$PERL" "$BOWTIE2" \
            -x "$SPIKE_INDEX" \
            -1 "$R1_PATH" \
            -2 "$R2_PATH" \
            -S "$spike_paired_sam" \
            --un-conc "${spike_unconc_pref}.fastq"


        echo "-> Aligning R1 singletons to spike RNA"


        if [ -s "$R1_UNPAIRED" ]; then

            "$PERL" "$BOWTIE2" \
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

            "$PERL" "$BOWTIE2" \
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
    # Align spike-unmapped reads directly to E. coli genome
    # ========================================================

    ecoli_paired_sam="${sample_dir}/${base}_eColi_paired.sam"

    ecoli_R1_unp_sam="${sample_dir}/${base}_eColi_R1_unpaired.sam"

    ecoli_R2_unp_sam="${sample_dir}/${base}_eColi_R2_unpaired.sam"


    spike_paired1="${spike_unconc_pref}.1.fastq"

    spike_paired2="${spike_unconc_pref}.2.fastq"


    if [[ -s "$ecoli_paired_sam" &&
          -f "$ecoli_R1_unp_sam" &&
          -f "$ecoli_R2_unp_sam" ]]; then

        echo "E. coli alignment already complete. Skipping."

    else

        # ----------------------------------------------------
        # Paired
        # ----------------------------------------------------

        echo "-> Aligning spike-unaligned paired reads to E. coli"


        if [[ -s "$spike_paired1" &&
              -s "$spike_paired2" ]]; then

            "$PERL" "$BOWTIE2" \
                -x "$ECOLI_INDEX" \
                -1 "$spike_paired1" \
                -2 "$spike_paired2" \
                -S "$ecoli_paired_sam" \
                --local

        else

            echo "No paired reads remain after spike alignment."

            : > "$ecoli_paired_sam"

        fi


        # ----------------------------------------------------
        # R1 singleton
        # ----------------------------------------------------

        echo "-> Aligning spike-unaligned R1 singletons to E. coli"


        if [ -s "$spike_R1_unp_fastq" ]; then

            "$PERL" "$BOWTIE2" \
                -x "$ECOLI_INDEX" \
                -U "$spike_R1_unp_fastq" \
                -S "$ecoli_R1_unp_sam" \
                --local

        else

            : > "$ecoli_R1_unp_sam"

        fi


        # ----------------------------------------------------
        # R2 singleton
        # ----------------------------------------------------

        echo "-> Aligning spike-unaligned R2 singletons to E. coli"


        if [ -s "$spike_R2_unp_fastq" ]; then

            "$PERL" "$BOWTIE2" \
                -x "$ECOLI_INDEX" \
                -U "$spike_R2_unp_fastq" \
                -S "$ecoli_R2_unp_sam" \
                --local

        else

            : > "$ecoli_R2_unp_sam"

        fi

    fi


    # --------------------------------------------------------
    # Sample finished successfully
    # --------------------------------------------------------

    echo
    echo "=== Finished $base ==="


    cleanup_lock


    echo "Lock released."


done


echo
echo "========================================"
echo "Step 04 completed."
echo "========================================"

echo "Alignment order:"
echo "  Spike -> E. coli"

echo

echo "Output directory:"
echo "$OUTPUT_DIR"

echo "========================================"
