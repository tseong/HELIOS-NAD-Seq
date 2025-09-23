#!/usr/bin/env python3
import glob
import os

# === CONFIGURATION ===
SAM_FOLDER = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/sam"
PATTERN    = "*eColi*eColi*.sam"

# --------------------------------------------------------------------------
def filter_sam_by_first_base(sam_path, out_path):
    """
    Copy headers and only those alignments where the first base of the read is A.
    This corresponds to the first base of the SEQ field in the SAM file.
    """
    with open(sam_path) as sam, open(out_path, 'w') as out:
        for line in sam:
            if line.startswith('@'):
                out.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 11:
                continue  # skip malformed lines

            seq = cols[9]  # 10th column = read sequence
            if seq.startswith('A'):
                out.write(line)

# --------------------------------------------------------------------------
def main():
    pattern   = os.path.join(SAM_FOLDER, PATTERN)
    sam_files = sorted(glob.glob(pattern))
    print(f"[1/2] Found {len(sam_files)} SAM files matching '{PATTERN}'")

    print("[2/2] Filtering reads starting with 'A'...")
    for sam in sam_files:
        base = os.path.basename(sam)
        out  = os.path.join(SAM_FOLDER, base.replace('.sam', '.Astart.sam'))
        print(f"    {base} ?~F~R {os.path.basename(out)}")
        filter_sam_by_first_base(sam, out)

    print("Done. Filtered SAMs are in:", SAM_FOLDER)

if __name__ == '__main__':
    main()
