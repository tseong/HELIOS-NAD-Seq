#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=90000
#SBATCH -t 8:00:00
#SBATCH -p cpu-single

set -euo pipefail

module load bio/bowtie2/2.4.5

# --- PATHS (edit indices as needed) ---
TRIM_DIR="/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/trimmed_trimmomatic"
SPIKE_INDEX="/home/hd/hd_hd/hd_uv268/software/bowtie2Index/spikeRna/spikeRnas"
TRNA_RRNA_INDEX="/home/hd/hd_hd/hd_uv268/software/bowtie2Index/trna_rrna/trna_rrna"   # <-- assume exists
ECOLI_INDEX="/home/hd/hd_hd/hd_uv268/software/bowtie2Index/eColiIndex/NC_00913.3_pUC19c"

cd "$TRIM_DIR"

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

    # ------------------------------------------------------------
    # 1) Align to SPIKE RNA (global), keep unaligned for next step
    # ------------------------------------------------------------
    spike_paired_sam="${base}_spike_paired.sam"
    spike_unpaired_r1_sam="${base}_spike_unpaired_r1.sam"
    spike_unpaired_r2_sam="${base}_spike_unpaired_r2.sam"

    spike_unconc_pref="${base}_unaligned_spike_paired"
    spike_R1_unp_fastq="${base}_unaligned_spike_R1_unpaired.fastq"
    spike_R2_unp_fastq="${base}_unaligned_spike_R2_unpaired.fastq"

    if [[ -f "$spike_paired_sam" && -f "$spike_unpaired_r1_sam" && -f "$spike_unpaired_r2_sam" ]]; then
        echo "-> Spike alignment exists, skipping"
    else
        echo "-> Aligning to Spike index (global mode)"
        bowtie2 \
            -x "$SPIKE_INDEX" \
            -1 "$r1" -2 "$r2" \
            -S "$spike_paired_sam" \
            --un-conc "${spike_unconc_pref}.fastq"

        bowtie2 \
            -x "$SPIKE_INDEX" \
            -U "$r1_unpaired" \
            -S "$spike_unpaired_r1_sam" \
            --un "$spike_R1_unp_fastq"

        bowtie2 \
            -x "$SPIKE_INDEX" \
            -U "$r2_unpaired" \
            -S "$spike_unpaired_r2_sam" \
            --un "$spike_R2_unp_fastq"
    fi

    # -------------------------------------------------------------------
    # 2) Align Spike-unaligned reads to tRNA + rRNA reference (depletion)
    #     Inputs:
    #       - Paired:  ${spike_unconc_pref}.1.fastq / .2.fastq
    #       - Singles: $spike_R1_unp_fastq / $spike_R2_unp_fastq
    #     Outputs (for next step):
    #       - Paired unaligned: ${rrna_unconc_pref}.1.fastq / .2.fastq
    #       - Singles unaligned: $rrna_R1_unp_fastq / $rrna_R2_unp_fastq
    # -------------------------------------------------------------------
    rrna_paired_sam="${base}_rrna_paired.sam"
    rrna_unpaired_r1_sam="${base}_rrna_unpaired_r1.sam"
    rrna_unpaired_r2_sam="${base}_rrna_unpaired_r2.sam"

    rrna_unconc_pref="${base}_unaligned_rrna_paired"
    rrna_R1_unp_fastq="${base}_unaligned_rrna_R1_unpaired.fastq"
    rrna_R2_unp_fastq="${base}_unaligned_rrna_R2_unpaired.fastq"

    spike_paired1="${spike_unconc_pref}.1.fastq"
    spike_paired2="${spike_unconc_pref}.2.fastq"

    if [[ -f "$rrna_paired_sam" && -f "$rrna_unpaired_r1_sam" && -f "$rrna_unpaired_r2_sam" ]]; then
        echo "-> tRNA+rRNA alignment exists, skipping"
    else
        echo "-> Aligning Spike-unaligned reads to tRNA+rRNA (depletion)"
        # Paired reads (from spike-unmapped)
        if [[ -s "$spike_paired1" && -s "$spike_paired2" ]]; then
            bowtie2 \
                -x "$TRNA_RRNA_INDEX" \
                -1 "$spike_paired1" -2 "$spike_paired2" \
                -S "$rrna_paired_sam" \
                --un-conc "${rrna_unconc_pref}.fastq"
        else
            # ensure empty placeholders exist for downstream step
            : > "${rrna_unconc_pref}.1.fastq"
            : > "${rrna_unconc_pref}.2.fastq"
        fi

        # Singletons
        if [[ -s "$spike_R1_unp_fastq" ]]; then
            bowtie2 \
                -x "$TRNA_RRNA_INDEX" \
                -U "$spike_R1_unp_fastq" \
                -S "$rrna_unpaired_r1_sam" \
                --un "$rrna_R1_unp_fastq"
        else
            : > "$rrna_R1_unp_fastq"
        fi

        if [[ -s "$spike_R2_unp_fastq" ]]; then
            bowtie2 \
                -x "$TRNA_RRNA_INDEX" \
                -U "$spike_R2_unp_fastq" \
                -S "$rrna_unpaired_r2_sam" \
                --un "$rrna_R2_unp_fastq"
        else
            : > "$rrna_R2_unp_fastq"
        fi
    fi

    # -------------------------------------------------------------------
    # 3) Align remaining unaligned reads to E. coli genome (local mode)
    #     Inputs come from the tRNA+rRNA UNALIGNED outputs
    # -------------------------------------------------------------------
    ecoli_paired_sam="${base}_eColi_paired.sam"
    ecoli_R1_unp_sam="${base}_eColi_R1_unpaired.sam"
    ecoli_R2_unp_sam="${base}_eColi_R2_unpaired.sam"

    ecoli_paired1="${rrna_unconc_pref}.1.fastq"
    ecoli_paired2="${rrna_unconc_pref}.2.fastq"

    if [[ -f "$ecoli_paired_sam" && -f "$ecoli_R1_unp_sam" && -f "$ecoli_R2_unp_sam" ]]; then
        echo "-> E. coli alignment exists, skipping"
    else
        echo "-> Aligning paired unaligned-from-rRNA/tRNA to E. coli (local mode)"
        if [[ -s "$ecoli_paired1" && -s "$ecoli_paired2" ]]; then
            bowtie2 \
                -x "$ECOLI_INDEX" \
                -1 "$ecoli_paired1" -2 "$ecoli_paired2" \
                -S "$ecoli_paired_sam" \
                --local
        else
            echo "   (no paired reads left after rRNA/tRNA depletion)"
            : > "$ecoli_paired_sam"
        fi

        echo "-> Aligning R1 singletons to E. coli"
        if [[ -s "$rrna_R1_unp_fastq" ]]; then
            bowtie2 \
                -x "$ECOLI_INDEX" \
                -U "$rrna_R1_unp_fastq" \
                -S "$ecoli_R1_unp_sam" \
                --local
        else
            : > "$ecoli_R1_unp_sam"
        fi

        echo "-> Aligning R2 singletons to E. coli"
        if [[ -s "$rrna_R2_unp_fastq" ]]; then
            bowtie2 \
                -x "$ECOLI_INDEX" \
                -U "$rrna_R2_unp_fastq" \
                -S "$ecoli_R2_unp_sam" \
                --local
        else
            : > "$ecoli_R2_unp_sam"
        fi
    fi

    echo "=== Done $base ==="
    echo
done
