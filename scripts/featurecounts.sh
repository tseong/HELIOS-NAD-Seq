#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

cd /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/sam_Astart

for file in *eColi*eColi*sam; do
    output_file="${file}.table"

    if [ ! -f "$output_file" ]; then
        echo "Processing $file..."

        # Check if the filename contains "_paired"
        if [[ "$file" == *_paired* ]]; then
            echo "Detected paired-end reads in $file. Using -p option."
            /home/hd/hd_hd/hd_uv268/software/subread-2.0.6-Linux-x86_64/bin/featureCounts -p \
                -a /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/data/gtf/GCF_000005845.2_ASM584v2_genomic_rna1_intergenic.gtf \
                -t intergenic \
                -o "$output_file" \
                "$file"
        else
            echo "Detected single-end reads in $file. Running without -p."
            /home/hd/hd_hd/hd_uv268/software/subread-2.0.6-Linux-x86_64/bin/featureCounts \
                -a /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/data/gtf/GCF_000005845.2_ASM584v2_genomic_rna1_intergenic.gtf \
                -t intergenic \
                -o "$output_file" \
                "$file"
        fi
    else
        echo "Skipping $file: $output_file already exists."
    fi
done
