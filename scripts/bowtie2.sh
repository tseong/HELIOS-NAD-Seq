#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

module load bio/bowtie2/2.4.5

cd /gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/trimmed_trimmomatic

for r1 in *R1_trimmed_paired.fastq; do
    echo "=== Sample: $r1 ==="

    # paired inputs
    r2="${r1//_R1_/_R2_}"
    if [[ ! -f "$r2" ]]; then
        echo "ERROR: Missing R2 for $r1" >&2
        continue
    fi

    # singleton inputs from trimming
    r1_unpaired="${r1/_paired.fastq/_unpaired.fastq}"
    r2_unpaired="${r1_unpaired//_R1_/_R2_}"

    base="${r1%_R1_trimmed_paired.fastq}"

    # 1) Spike RNA alignment
    spike_paired_sam="${base}_spike_paired.sam"
    spike_unpaired_r1_sam="${base}_spike_unpaired_r1.sam"
    spike_unpaired_r2_sam="${base}_spike_unpaired_r2.sam"
    spike_unconc_pref="${base}_unaligned_spike_paired"
    spike_R1_unp_fastq="${base}_unaligned_spike_R1_unpaired.fastq"
    spike_R2_unp_fastq="${base}_unaligned_spike_R2_unpaired.fastq"

    # Skip if already done
    if [[ -f "$spike_paired_sam" || -f "$spike_unpaired_r1_sam" || -f "$spike_unpaired_r2_sam" ]]; then
        echo "-> Spike alignment exists, skipping"
       continue
    fi

    echo "-> Aligning to Spike index (global mode)"
    bowtie2 \
        -x /home/hd/hd_hd/hd_uv268/software/bowtie2Index/spikeRna/spikeRnas \
        -1 "$r1" -2 "$r2" \
        -S "$spike_paired_sam" \
        --un-conc "${spike_unconc_pref}.fastq"

    bowtie2 -x /home/hd/hd_hd/hd_uv268/software/bowtie2Index/spikeRna/spikeRnas \
        -U $r1_unpaired \
        -S $spike_unpaired_r1_sam \
        --un $spike_R1_unp_fastq

     Align unpaired R2 reads
    bowtie2 -x /home/hd/hd_hd/hd_uv268/software/bowtie2Index/spikeRna/spikeRnas \
        -U $r2_unpaired \
        -S $spike_unpaired_r2_sam \
        --un $spike_R2_unp_fastq       
    # 2) E coli alignment of the reads that failed to map to Spike

    ecoli_paired_sam="${base}_eColi_paired.sam"
    ecoli_R1_unp_sam="${base}_eColi_R1_unpaired.sam"
    ecoli_R2_unp_sam="${base}_eColi_R2_unpaired.sam"
    ecoli_paired1="${spike_unconc_pref}.1.fastq"
    ecoli_paired2="${spike_unconc_pref}.2.fastq"

    if [[ -f "$ecoli_paired_sam" || -f "$ecoli_R1_unp_sam" || -f "$ecoli_R2_unp_sam" ]]; then
        echo "-> Spike alignment exists, skipping"
       continue
    fi

    echo "-> Aligning paired unaligned?~@~Pfrom?~@~PSpike to E?| coli (local mode)"
    bowtie2 \
        -x /home/hd/hd_hd/hd_uv268/software/bowtie2Index/eColiIndex/NC_00913.3_pUC19c \
        -1 "$ecoli_paired1" -2 "$ecoli_paired2" \
        -S "$ecoli_paired_sam" \
        --local

    echo "-> Aligning R1 singletons to E?| coli"
    bowtie2 \
        -x /home/hd/hd_hd/hd_uv268/software/bowtie2Index/eColiIndex/NC_00913.3_pUC19c \
        -U "$spike_R1_unp_fastq" \
        -S "$ecoli_R1_unp_sam" \
        --local

    echo "-> Aligning R2 singletons to E?| coli"
    bowtie2 \
        -x /home/hd/hd_hd/hd_uv268/software/bowtie2Index/eColiIndex/NC_00913.3_pUC19c \
        -U "$spike_R2_unp_fastq" \
        -S "$ecoli_R2_unp_sam" \
        --local

    echo "=== Done $base ==="
    echo
done
