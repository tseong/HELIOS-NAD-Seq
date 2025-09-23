#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

# Trimming with Trimmomatic
faPara=('1:30:15')

end_dir=/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/cutadapt
cd $end_dir

for file in *R1_001_trimmed.fastq; do
    # Define R1 and R2 filenames based on the R1 file name
    R1_file="$file"
    R2_file="${file/_R1_/_R2_}"

    # Define output file names for paired and unpaired reads
    base_name="${file%_R1_001_trimmed.fastq}"
    R1_paired="${base_name}_R1_trimmed_paired.fastq"
    R1_unpaired="${base_name}_R1_trimmed_unpaired.fastq"
    R2_paired="${base_name}_R2_trimmed_paired.fastq"
    R2_unpaired="${base_name}_R2_trimmed_unpaired.fastq"

    if [[ -f "$R1_paired" && -f "$R1_unpaired" && -f "$R2_paired" && -f "$R2_unpaired" ]]; then
        echo "Trimmomatic output for $base_name already exists. Skipping..."
        continue
    fi

    # Run Trimmomatic in PE mode
    /home/hd/hd_hd/hd_uv268/software/jdk1.8.0_361/bin/java -jar /home/hd/hd_hd/hd_uv268/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
        -trimlog /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/logs/200425_trimmomatic.trimlog \
        "$R1_file" "$R2_file" \
        "$R1_paired" "$R1_unpaired" "$R2_paired" "$R2_unpaired" \
        ILLUMINACLIP:/home/hd/hd_hd/hd_uv268/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:$faPara SLIDINGWINDOW:4:20 MINLEN:18

done
