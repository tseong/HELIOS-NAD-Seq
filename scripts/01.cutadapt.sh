#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

# ============================================================
# Usage:
#   bash 01.cutadapt.sh <FASTQ_DIR> [OUTPUT_DIR]
#
# Example:
#   bash 01.cutadapt.sh /path/to/fastq
#
#   or:
#   bash 01.cutadapt.sh /path/to/fastq /path/to/cutadapt
#
# If OUTPUT_DIR is not provided, a directory named "cutadapt"
# will be created next to the FASTQ directory.
# ============================================================


# Check that an input directory was provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <FASTQ_DIR> [OUTPUT_DIR]"
    exit 1
fi


# Input FASTQ directory
start_dir="$1"

# Remove trailing slash, if present
start_dir="${start_dir%/}"


# Check that input directory exists
if [ ! -d "$start_dir" ]; then
    echo "ERROR: FASTQ directory does not exist:"
    echo "$start_dir"
    exit 1
fi


# Set output directory
if [ $# -ge 2 ]; then
    end_dir="$2"
else
    parent_dir="$(dirname "$start_dir")"
    end_dir="${parent_dir}/cutadapt"
fi


# Create output directory if necessary
mkdir -p "$end_dir"


echo "FASTQ input directory: $start_dir"
echo "Cutadapt output directory: $end_dir"


# Move to FASTQ directory
cd "$start_dir" || exit 1


# Find paired-end FASTQ files
forward_reads=(*R1*.fastq)
reverse_reads=(*R2*.fastq)


# Check that FASTQ files were found
if [ ! -e "${forward_reads[0]}" ]; then
    echo "ERROR: No R1 FASTQ files found in $start_dir"
    exit 1
fi

if [ ! -e "${reverse_reads[0]}" ]; then
    echo "ERROR: No R2 FASTQ files found in $start_dir"
    exit 1
fi


# Check that the number of R1 and R2 files matches
if [ "${#forward_reads[@]}" -ne "${#reverse_reads[@]}" ]; then
    echo "ERROR: Number of R1 and R2 FASTQ files does not match."
    echo "R1 files: ${#forward_reads[@]}"
    echo "R2 files: ${#reverse_reads[@]}"
    exit 1
fi

i=0

for forward_read in "${forward_reads[@]}"; do

    # Define paired reads
    paired_read=("${forward_reads[$i]}" "${reverse_reads[$i]}")

    R1_base="$(basename "${paired_read[0]}")"
    R2_base="$(basename "${paired_read[1]}")"

    echo "========================================"
    echo "Checking pair:"
    echo "  R1: ${paired_read[0]}"
    echo "  R2: ${paired_read[1]}"

    # --------------------------------------------------------
    # Check whether all expected barcode outputs already exist
    # --------------------------------------------------------

    all_outputs_exist=true

    for bc in bc01 bc02 bc03 bc04 bc05 bc06 bc07 bc08; do

        # Note: R2 is primary input (-o), R1 is paired input (-p)
        output_R2="$end_dir/${bc}_${R2_base}"
        output_R1="$end_dir/${bc}_${R1_base}"

        if [[ ! -s "$output_R1" || ! -s "$output_R2" ]]; then
            all_outputs_exist=false
            break
        fi

    done

    if [[ "$all_outputs_exist" == true ]]; then

        echo "All bc01-bc08 output files already exist."
        echo "Skipping this FASTQ pair."

        i=$((i+1))
        continue

    fi

    echo "Output incomplete or not found."
    echo "Running Cutadapt..."

    cutadapt \
        -j 0 \
        -O 12 \
        -g bc01=TCAAGTNNNNNNG \
        -g bc01=TCAAGTNNNNNNNG \
        -g bc02=CAGCGTNNNNNNG \
        -g bc02=CAGCGTNNNNNNNG \
        -g bc03=ACCGGTNNNNNNG \
        -g bc03=ACCGGTNNNNNNNG \
        -g bc04=ATGAGTNNNNNNG \
        -g bc04=ATGAGTNNNNNNNG \
        -g bc05=GTTCGTNNNNNNG \
        -g bc05=GTTCGTNNNNNNNG \
        -g bc06=TGCTGTNNNNNNG \
        -g bc06=TGCTGTNNNNNNNG \
        -g bc07=TATGGTNNNNNNG \
        -g bc07=TATGGTNNNNNNNG \
        -g bc08=CTATGTNNNNNNG \
        -g bc08=CTATGTNNNNNNNG \
        -A CNNNNNNACTTGA \
        -A CNNNNNNACGCTG \
        -A CNNNNNNACCGGT \
        -A CNNNNNNACTCAT \
        -A CNNNNNNACGAAC \
        -A CNNNNNNACAGCA \
        -A CNNNNNNACCATA \
        -A CNNNNNNACATAG \
        -A CNNNNNNNACTTGA \
        -A CNNNNNNNACGCTG \
        -A CNNNNNNNACCGGT \
        -A CNNNNNNNACTCAT \
        -A CNNNNNNNACGAAC \
        -A CNNNNNNNACAGCA \
        -A CNNNNNNNACCATA \
        -A CNNNNNNNACATAG \
        -o "$end_dir/{name}_${R2_base}" \
        -p "$end_dir/{name}_${R1_base}" \
        "${paired_read[1]}" "${paired_read[0]}"

    i=$((i+1))

done

echo "Demultiplexing complete."
echo "Output directory: $end_dir"
