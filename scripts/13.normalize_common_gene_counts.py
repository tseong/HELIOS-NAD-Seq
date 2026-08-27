#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 13: Normalize common variable-3' gene counts
#
# This step normalizes common-gene total counts using the
# assigned-read depth at each time point.
#
# Normalization factor:
#
#     tp1 assigned reads / tpX assigned reads
#
#
# Expected input structure:
#
# count_tables/
# └── variable_3prime/
#     ├── assigned_reads.csv
#     │
#     └── common_genes_tp_least10_assigned_read_normalized/
#         └── common_genes_total_counts_summary.csv
#
#
# Output:
#
# common_genes_tp_least10_assigned_read_normalized/
# └── common_genes_total_counts_summary_assigned_reads_normalized.csv
#
#
# Usage:
#
#   python scripts/13.normalize_common_gene_counts.py \
#       <VARIABLE_3PRIME_DIR>
#
#
# Example:
#
#   python scripts/13.normalize_common_gene_counts.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables/variable_3prime
#
# ============================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 13: normalize common variable-3' "
            "gene counts using assigned-read depth and calculate "
            "the percent coefficient of variation across time points."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Variable 3' count-table directory containing "
            "assigned_reads.csv and the common-gene normalization "
            "subdirectory."
        )
    )


    parser.add_argument(
        "--start-tp",
        type=int,
        default=1,
        help="First time point to normalize. Default: 1."
    )


    parser.add_argument(
        "--end-tp",
        type=int,
        default=16,
        help="Last time point to normalize. Default: 16."
    )


    parser.add_argument(
        "--reference-tp",
        default="tp1",
        help=(
            "Reference time point used for normalization. "
            "Default: tp1."
        )
    )


    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help=(
            "Minimum-count threshold used in the previous preparation "
            "step. This determines the expected common-gene folder name. "
            "Default: 10."
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
            f"Input directory does not exist:\n"
            f"{input_dir}"
        )


    # --------------------------------------------------------
    # Validate time points
    # --------------------------------------------------------

    if args.start_tp < 1:

        raise ValueError(
            "--start-tp must be >= 1."
        )


    if args.end_tp < args.start_tp:

        raise ValueError(
            "--end-tp must be >= --start-tp."
        )


    tp_columns = [
        f"tp{i}"
        for i in range(
            args.start_tp,
            args.end_tp + 1
        )
    ]


    if args.reference_tp not in tp_columns:

        raise ValueError(
            f"--reference-tp must be one of: "
            f"{', '.join(tp_columns)}"
        )


    # --------------------------------------------------------
    # Expected input/output files
    # --------------------------------------------------------

    read_depth_file = (
        input_dir
        / "assigned_reads.csv"
    )


    normalization_dir = (
        input_dir
        / (
            f"common_genes_tp_least"
            f"{args.min_count}_assigned_read_normalized"
        )
    )


    counts_file = (
        normalization_dir
        / "common_genes_total_counts_summary.csv"
    )


    output_file = (
        normalization_dir
        / "common_genes_total_counts_summary_assigned_reads_normalized.csv"
    )


    # --------------------------------------------------------
    # Validate input files
    # --------------------------------------------------------

    if not read_depth_file.is_file():

        raise FileNotFoundError(
            "Assigned-read file does not exist:\n"
            f"{read_depth_file}"
        )


    if not counts_file.is_file():

        raise FileNotFoundError(
            "Common-gene count file does not exist:\n"
            f"{counts_file}"
        )


    normalization_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 13 - Assigned-read normalization")
    print("========================================")

    print(
        f"Variable 3' directory: {input_dir}"
    )

    print(
        f"Count file:            {counts_file}"
    )

    print(
        f"Assigned-read file:    {read_depth_file}"
    )

    print(
        f"Time points:           "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        f"Reference time point:  {args.reference_tp}"
    )

    print(
        f"Output:                {output_file}"
    )

    print("========================================")


    # ========================================================
    # Load input data
    # ========================================================

    counts_df = pd.read_csv(
        counts_file
    )


    read_depth_df = pd.read_csv(
        read_depth_file
    )


    # --------------------------------------------------------
    # Validate time-point columns
    # --------------------------------------------------------

    missing_count_tp = [
        tp
        for tp in tp_columns
        if tp not in counts_df.columns
    ]


    if missing_count_tp:

        raise ValueError(
            "Missing time-point columns in count file: "
            + ", ".join(missing_count_tp)
        )


    missing_depth_tp = [
        tp
        for tp in tp_columns
        if tp not in read_depth_df.columns
    ]


    if missing_depth_tp:

        raise ValueError(
            "Missing time-point columns in assigned_reads.csv: "
            + ", ".join(missing_depth_tp)
        )


    if len(read_depth_df) < 1:

        raise ValueError(
            "assigned_reads.csv contains no data rows."
        )


    # --------------------------------------------------------
    # Convert values to numeric
    # --------------------------------------------------------

    counts_df[
        tp_columns
    ] = counts_df[
        tp_columns
    ].apply(
        pd.to_numeric,
        errors="raise"
    )


    read_depths = pd.to_numeric(
        read_depth_df.loc[
            0,
            tp_columns
        ],
        errors="raise"
    )


    # --------------------------------------------------------
    # Validate assigned-read depths
    # --------------------------------------------------------

    if (
        read_depths <= 0
    ).any():

        invalid_tp = list(
            read_depths[
                read_depths <= 0
            ].index
        )

        raise ValueError(
            "Assigned-read depth must be > 0 for all time points. "
            "Invalid values detected for: "
            + ", ".join(invalid_tp)
        )


    # ========================================================
    # Calculate normalization factors
    # ========================================================

    reference_depth = (
        read_depths[
            args.reference_tp
        ]
    )


    norm_factors = (
        reference_depth
        / read_depths
    )


    print()
    print("Normalization factors:")


    for tp in tp_columns:

        print(
            f"  {tp}: "
            f"{norm_factors[tp]:.6f}"
        )


    # ========================================================
    # Normalize counts
    # ========================================================

    normalized_df = (
        counts_df.copy()
    )


    for tp in tp_columns:

        normalized_df[
            tp
        ] = (
            normalized_df[
                tp
            ]
            * norm_factors[
                tp
            ]
        )


    # ========================================================
    # Calculate %CV across time points
    #
    # %CV = standard deviation / mean * 100
    # ========================================================

    means = normalized_df[
        tp_columns
    ].mean(
        axis=1
    )


    stds = normalized_df[
        tp_columns
    ].std(
        axis=1
    )


    normalized_df[
        "%CV"
    ] = np.where(
        means != 0,
        (
            stds
            / means
        ) * 100,
        np.nan
    )


    # --------------------------------------------------------
    # Save output
    # --------------------------------------------------------

    normalized_df.to_csv(
        output_file,
        index=False
    )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("========================================")
    print("Step 13 completed.")
    print("----------------------------------------")

    print(
        f"Genes normalized: "
        f"{len(normalized_df)}"
    )

    print(
        f"Reference depth "
        f"({args.reference_tp}): "
        f"{reference_depth}"
    )

    print()

    print(
        "Normalized output:"
    )

    print(
        output_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
