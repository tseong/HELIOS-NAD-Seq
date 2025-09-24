#!/usr/bin/env python3
import subprocess
from pathlib import Path
from collections import defaultdict

# === CONFIG ===
INPUT_DIR = Path("/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/sam/intermediate")
LOGO_DIR  = Path("/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/sam/270525_logos_by_rel_position_0_final")
LOGO_DIR.mkdir(parents=True, exist_ok=True)

LOGO_FORMAT = "pdf"   # vector output: pdf/png/svg/eps
DPI = 600

# Expected window on the sequences already present in the FASTA:
#   -40 .. 0 .. +4  (inclusive)  => length 45
WINDOW_START = -40
WINDOW_END   = 4
EXPECTED_LEN = WINDOW_END - WINDOW_START + 1  # 45

# Build custom labels for WebLogo: -40..-1, +1..+5 (0 -> +1)
labels = [str(i) for i in range(WINDOW_START, 0)] + [f"+{i}" for i in range(1, WINDOW_END + 2)]
ANNOT_STR = ",".join(labels)  # must have length == EXPECTED_LEN (45)

def make_weighted_fasta(sample: str, chrom: str, seqs: list[str], out_path: Path):
    with open(out_path, "w") as fout:
        for i, seq in enumerate(seqs):
            fout.write(f">{sample}_{chrom}_{i}\n{seq}\n")

# Process each FASTA file independently
for fasta_path in INPUT_DIR.glob("*.fa"):
    sample = fasta_path.stem  # e.g., 'top10_relpos_bc01_tp1'
    print(f"Processing {fasta_path.name}")

    # 1) Read all sequences into per-chrom bins, weighted by count
    chrom_seqs = defaultdict(list)
    with open(fasta_path) as fin:
        current_count = 1
        current_chrom = None
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                # header format expected: COUNT::CHROM:POS-POS(STRAND)
                parts = header.split("::", 1)
                try:
                    current_count = int(parts[0])
                except Exception:
                    current_count = 1
                rest = parts[1] if len(parts) > 1 else ""
                current_chrom = rest.split(":", 1)[0] if ":" in rest else "unknown"
            else:
                seq = line.strip().upper()
                # replicate the sequence by its count
                chrom_seqs[current_chrom].extend([seq] * current_count)
                  # 2) For each chromosome present, write a weighted FASTA & generate a logo
    for chrom, seqs in chrom_seqs.items():
        if not seqs:
            continue

        # Sanity: check sequence length
        seq_lens = {len(s) for s in seqs}
        use_annotate = True
        if seq_lens != {EXPECTED_LEN}:
            use_annotate = False
            print(f"  [WARN] {sample}|{chrom}: sequence lengths {sorted(seq_lens)} "
                  f"!= expected {EXPECTED_LEN}. Using default indexing instead of custom labels.")

        # write weighted FASTA
        weighted_fa = LOGO_DIR / f"{sample}_{chrom}_weighted.fa"
        make_weighted_fasta(sample, chrom, seqs, weighted_fa)

        # generate logo
        out_logo = LOGO_DIR / f"{sample}_{chrom}.{LOGO_FORMAT}"
        cmd = [
            "weblogo",
            "-f", str(weighted_fa),
            "-o", str(out_logo),
            "--format", LOGO_FORMAT,
            "--resolution", str(DPI),
            "--title", f"{sample} | {chrom}  (TSS labeled +1)",
            "--xlabel", "Position relative to TSS",
            "--units", "bits",
            "--size", "large",
            "--number-interval", "5",
        ]

        # Use custom labels (keeps negatives as-is; maps 0->+1, 1->+2, ..., 4->+5)
        if use_annotate:
            cmd += ["--annotate", ANNOT_STR]
        else:
            # Fall back to natural indexing starting at WINDOW_START
            # (This won?~@~Yt re-label 0?~F~R+1; it just avoids mismatched annotate length.)
            cmd += ["--first-index", str(WINDOW_START)]

        try:
            subprocess.run(cmd, check=True)
            print(f"  Created logo for {chrom}: {out_logo.name}")
        except subprocess.CalledProcessError as e:
            print(f"  Error generating logo for {chrom}: {e}")

        # clean up weighted FASTA
        weighted_fa.unlink()

print("All logos generated. Labels use +1 at TSS (0 shifted) when sequence length matches -40..+4.")
