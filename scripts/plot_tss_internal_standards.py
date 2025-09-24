#!/usr/bin/env python3
import os
import glob
import argparse
import matplotlib.pyplot as plt
from collections import Counter

def build_tss_counts(sam_files):
    """
    Build per-reference TSS counts (position of alignment start) from SAM files.
    Returns a dict: { rname: Counter({pos: count, ...}) }.
    """
    tss = {}
    for path in sam_files:
        with open(path, "r") as f:
            for line in f:
                if line.startswith("@"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 4:
                    continue
                rname = fields[2]
                pos   = int(fields[3])  # 1-based
                cnt = tss.setdefault(rname, Counter())
                cnt[pos] += 1
    return tss

def plot_tss(tss_counts, outdir, bc, tp):
    """
    Plot bar chart of TSS counts for each reference in tss_counts.
    Save as PNG if majority at position 1, else PDF.
    """
    os.makedirs(outdir, exist_ok=True)
    for rname, counts in tss_counts.items():
        total = sum(counts.values())
        first = counts.get(1, 0)
        majority = (first > total / 5)
        if majority:
            fname = f"{rname}_tss_majority.png"
        else:
            fname = f"{rname}_tss.png"
        outpath = os.path.join(outdir, fname)

        positions = sorted(counts)
        values = [counts[p] for p in positions]

        plt.figure(figsize=(8, 3))
        plt.bar(positions, values)
        plt.title(f"TSS distribution for {rname} ({bc} {tp})")
        plt.xlabel("TSS position (1-based)")
        plt.ylabel("Read count")
        plt.tight_layout()
        plt.savefig(outpath)
        plt.close()
      def main():
    parser = argparse.ArgumentParser(
        description="Plot TSS distributions for each barcode (bc01?~@~Sbc08) and "
                    "timepoint (tp1?~@~Stp16) from SAM files.")
    parser.add_argument("sam_dir",
                        help="Directory containing SAM files and where output subfolders will be created")
    args = parser.parse_args()

    sam_dir = args.sam_dir
    barcodes = [f"bc{i:02d}" for i in range(1, 9)]
    timepoints = [f"tp{i}" for i in range(1, 17)]

    for tp in timepoints:
        for bc in barcodes:
            pattern = os.path.join(sam_dir, f"{bc}*{tp}_*spike*.sam")
            sam_files = glob.glob(pattern)
            if not sam_files:
                continue

            print(f"Processing {bc} {tp}: {len(sam_files)} files")
            tss_counts = build_tss_counts(sam_files)

            outdir = os.path.join(sam_dir, tp, bc)
            plot_tss(tss_counts, outdir, bc, tp)
            print(f"  -> Plots saved in {outdir}")

if __name__ == "__main__":
    main()
