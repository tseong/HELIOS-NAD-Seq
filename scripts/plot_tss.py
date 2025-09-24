import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import os

# Directory containing merged TSV files
INPUT_DIR = Path("/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/sam/intermediate")
OUTPUT_DIR = INPUT_DIR.parent / "tss_positions"
# Attempt to create output dir, ignore permissions if it fails
try:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
except PermissionError:
    pass

# pattern to extract barcode & timepoint
pattern = re.compile(r"read_starts_relative_position_bc(\d{2})_eColi_tp(\d{1,2})\.bed$")

# chromosomes to plot separately
chroms = ["NC_000913.3", "puc19C"]

for tsv in INPUT_DIR.glob("read_starts_relative_position_bc*_eColi_tp*.bed"):
    m = pattern.search(tsv.name)
    if not m:
        continue
    bc, tp = m.groups()
    prefix = f"bc{bc}_tp{tp}"
    print(f"Processing {tsv.name} ?~F~R {prefix}")

    # read files with columns: chrom, start, end, rel_pos, count, strand
    df = pd.read_csv(tsv, sep="\t", header=0,
                     usecols=["chrom","start","end","rel_pos","count","strand"])
    df["rel_pos"] = pd.to_numeric(df["rel_pos"], errors="raise")
    df["count"]   = pd.to_numeric(df["count"],   errors="raise")

    # for each chromosome, make its own PDF
    for chrom in chroms:
        sub = df[df.chrom == chrom]
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(10,4))

        # sort by rel_pos and plot
        sub = sub.sort_values("rel_pos")
        ax.bar(sub.rel_pos, sub["count"], width=0.8, label=chrom, color="#1F77B4")

        ax.set_xlabel("Position relative to TSS")
        ax.set_ylabel("Read count")
        ax.set_title(f"{chrom} read start distribution: {prefix}")

        # x-axis from -30 to +30
        ticks = list(range(-30, 31, 5))
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_xlim(-30, 30)
              # compute max only over the window rel_pos ?~H~H [?~@~S30, +30]
        window = sub[(sub.rel_pos >= -30) & (sub.rel_pos <= 30)]
        if not window.empty:
            max_count = window["count"].max()
        else:
            # fallback if nothing falls in the window
            max_count = sub["count"].max()

        ax.set_ylim(0, max_count * 1.1)

        ax.legend().set_visible(False)
        plt.tight_layout()

        out_pdf = OUTPUT_DIR / f"read_start_{prefix}_{chrom}.pdf"
        try:
            plt.savefig(out_pdf, format="pdf")
        except PermissionError:
            # fallback to current dir
            fallback = Path(f"read_start_{prefix}_{chrom}.pdf")
            plt.savefig(fallback, format="pdf")
        plt.close(fig)
