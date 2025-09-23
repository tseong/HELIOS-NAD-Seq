#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

cd /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/fastq

mkdir /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/cutadapt

start_dir=/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/fastq
end_dir=/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/cutadapt

forward_reads=(*R1*.fastq)
reverse_reads=(*R2*.fastq)

i=0

for forward_read in ${forward_reads[@]}; do
        # Define paired reads correctly
        paired_read=("${forward_reads[$i]}" "${reverse_reads[$i]}")

        # Define output filenames
        output_forward="$end_dir/bc01_$(basename "${paired_read[0]}")"
        output_reverse="$end_dir/bc01_$(basename "${paired_read[1]}")"

        # Check if output files exist
        if [[ -f "$output_forward" && -f "$output_reverse" ]]; then
                echo "Output files $output_forward and $output_reverse already exist. Skipping..."
        else
                echo "Processing ${paired_read[0]} and ${paired_read[1]} because output does not exist."
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
        -o "$end_dir/{name}_$(basename "${paired_read[1]}")" \
        -p "$end_dir/{name}_$(basename "${paired_read[0]}")" \
        "${paired_read[1]}" "${paired_read[0]}"
        #demultiplexing is done for the R2 file specified as an R1 file because cutadapt only seems to support demultiplexing R1 files

        fi

        i=$((i+1))

done
