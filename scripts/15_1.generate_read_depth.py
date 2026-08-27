#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 15_1: Generate total read depth by time point
#
# This script calculates total read depth for tp1-tp16 from
# FASTQ files and generates:
#
#   read_depth_by_timepoint_all_reads.csv
#
#
# Expected input:
#
# trimmomatic/
# ├── bc01_tp1_R1_trimmed_paired.fastq
# ├── bc01_tp1_R1_trimmed_unpaired.fastq
# ├── bc01_tp1_R2_trimmed_paired.fastq
# ├── bc01_tp1_R2_trimmed_unpaired.fastq
# ├── ...
#
#
# The script counts FASTQ records and sums them by time point.
#
# Paired reads are counted from R1 paired files only, so each
# fragment is counted once rather than twice.
#
# R1 and R2 unpaired reads are both counted because they
# represent separate surviving singleton reads.
#
#
# Output:
#
# trimmomatic/
# └── read_depth_by_timepoint_all_reads.csv
#
#
# Usage:
#
# python scripts/15_1.generate_read_depth.py \
#     <TRIMMOMATIC_DIR>
#
#
# Example:
#
# python scripts/15_1.generate_read_depth.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/trimmomatic
#
# ============================================================


def count_fastq_reads(fastq_path):
    """
    Count FASTQ records.

    FASTQ contains four lines per read.
    """

    line_count = 0

    with fastq_path.open("r") as handle:

        for _ in handle:
            line_count += 1

    if line_count % 4 != 0:

        raise ValueError(
            f"FASTQ file does not contain a multiple of 4 lines:\n"
            f"{fastq_path}\n"
            f"Lines: {line_count}"
        )

    return line_count // 4


def extract_timepoint(filename):
    """
    Extract tpX from filename.

    Example:
        bc01_tp10_R1_trimmed_paired.fastq
        -> tp10
    """

    match = re.search(
        r'(tp\d+)',
        filename
    )

    if match:
        return match.group(1)

    return None


def tp_sort_key(tp):
    """
    Sort tp1, tp2, ..., tp16 numerically.
    """

    return int(
        tp.replace("tp", "")
    )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 15_1: calculate total read "
            "depth per time point from Trimmomatic FASTQ files."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Trimmomatic output directory containing paired "
            "and unpaired FASTQ files."
        )
    )


    parser.add_argument(
        "--start-tp",
        type=int,
        default=1,
        help="First time point. Default: 1."
    )


    parser.add_argument(
        "--end-tp",
        type=int,
        default=16,
        help="Last time point. Default: 16."
    )


    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output CSV. Default: "
            "<input_dir>/read_depth_by_timepoint_all_reads.csv"
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Resolve paths
    # --------------------------------------------------------

    input_dir = Path(
        args.input_dir
    ).expanduser().resolve()


    if not input_dir.is_dir():

        raise FileNotFoundError(
            f"Trimmomatic directory does not exist:\n"
            f"{input_dir}"
        )


    if args.start_tp < 1:

        raise ValueError(
            "--start-tp must be >= 1."
        )


    if args.end_tp < args.start_tp:

        raise ValueError(
            "--end-tp must be >= --start-tp."
        )


    if args.output:

        output_file = Path(
            args.output
        ).expanduser().resolve()

    else:

        output_file = (
            input_dir
            / "read_depth_by_timepoint_all_reads.csv"
        )


    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Initialize time points
    # --------------------------------------------------------

    tp_names = [
        f"tp{i}"
        for i in range(
            args.start_tp,
            args.end_tp + 1
        )
    ]


    read_depths = {
        tp: 0
        for tp in tp_names
    }


    # --------------------------------------------------------
    # Find FASTQ files
    #
    # Count:
    #
    #   R1 paired
    #   R1 unpaired
    #   R2 unpaired
    #
    # Do NOT count R2 paired because paired fragments would
    # otherwise be counted twice.
    # --------------------------------------------------------

    paired_r1_files = sorted(
        input_dir.glob(
            "*R1_trimmed_paired.fastq"
        )
    )


    unpaired_r1_files = sorted(
        input_dir.glob(
            "*R1_trimmed_unpaired.fastq"
        )
    )


    unpaired_r2_files = sorted(
        input_dir.glob(
            "*R2_trimmed_unpaired.fastq"
        )
    )


    all_files = (
        paired_r1_files
        + unpaired_r1_files
        + unpaired_r2_files
    )


    if not all_files:

        raise FileNotFoundError(
            "No expected Trimmomatic FASTQ files were found in:\n"
            f"{input_dir}"
        )


    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 15_1 - Generate read depth")
    print("========================================")
    print(f"Input directory: {input_dir}")
    print(
        f"Time points:     "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )
    print(f"Output:          {output_file}")
    print("----------------------------------------")
    print(f"R1 paired files:   {len(paired_r1_files)}")
    print(f"R1 unpaired files: {len(unpaired_r1_files)}")
    print(f"R2 unpaired files: {len(unpaired_r2_files)}")
    print("========================================")


    # ========================================================
    # Count reads
    # ========================================================

    file_count = 0


    for fastq_path in all_files:

        tp = extract_timepoint(
            fastq_path.name
        )


        if tp is None:

            print(
                f"WARNING: Could not detect time point in "
                f"{fastq_path.name}. Skipping."
            )

            continue


        if tp not in read_depths:

            continue


        n_reads = count_fastq_reads(
            fastq_path
        )


        read_depths[
            tp
        ] += n_reads


        print(
            f"{fastq_path.name}: "
            f"{n_reads} reads -> {tp}"
        )


        file_count += 1


    # --------------------------------------------------------
    # Check all requested time points
    # --------------------------------------------------------

    missing_tp = [
        tp
        for tp in tp_names
        if read_depths[tp] == 0
    ]


    if missing_tp:

        print()
        print(
            "WARNING: Zero total reads detected for: "
            + ", ".join(
                missing_tp
            )
        )


    # ========================================================
    # Write one-row CSV
    # ========================================================

    output_df = pd.DataFrame(
        [
            {
                tp: read_depths[tp]
                for tp in tp_names
            }
        ]
    )


    output_df.to_csv(
        output_file,
        index=False
    )


    # ========================================================
    # Summary
    # ========================================================

    print()
    print("========================================")
    print("Step 15_1 completed.")
    print("----------------------------------------")
    print(f"FASTQ files counted: {file_count}")
    print()

    for tp in sorted(
        tp_names,
        key=tp_sort_key
    ):

        print(
            f"{tp}: {read_depths[tp]}"
        )

    print()
    print("Output:")
    print(
        output_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
