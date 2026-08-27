#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Prepare variable-3' normalization inputs
#
# Generates:
#
#   1. assigned_reads.csv
#
#   2. common_genes_total_counts_summary.csv
#
#
# Expected input:
#
# count_tables/
# └── variable_3prime/
#     ├── tp1/
#     │   └── merged_by_barcode_Astart_readCount.csv
#     ├── tp2/
#     │   └── merged_by_barcode_Astart_readCount.csv
#     └── ...
#
#
# NOTE:
# The filename may still contain "Astart" because it was
# inherited from the earlier pipeline naming, although the
# variable_3prime counts themselves came from the original
# Step 04 alignments.
#
#
# Output:
#
# variable_3prime/
# ├── assigned_reads.csv
# │
# └── common_genes_tp_least10_assigned_read_normalized/
#     └── common_genes_total_counts_summary.csv
#
#
# Usage:
#
# python scripts/11_2.prepare_variable3prime_normalization.py \
#     <VARIABLE_3PRIME_COUNT_DIR>
#
#
# Example:
#
# python scripts/prepare_variable3prime_normalization.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables/variable_3prime
#
# ============================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare assigned-read depths and common-gene "
            "total-count matrices for variable 3' normalization."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Variable 3' count-table directory containing "
            "tp1, tp2, ... subdirectories."
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
        "--min-count",
        type=int,
        default=10,
        help=(
            "Minimum total count required for a gene at every "
            "time point to be included in the common-gene set. "
            "Default: 10."
        )
    )


    parser.add_argument(
        "--input-name",
        default="merged_by_barcode_Astart_readCount.csv",
        help=(
            "Merged count filename in each tp directory. "
            "Default: merged_by_barcode_Astart_readCount.csv"
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


    if args.start_tp < 1:

        raise ValueError(
            "--start-tp must be >= 1."
        )


    if args.end_tp < args.start_tp:

        raise ValueError(
            "--end-tp must be >= --start-tp."
        )


    if args.min_count < 0:

        raise ValueError(
            "--min-count must be >= 0."
        )


    # --------------------------------------------------------
    # Output locations
    # --------------------------------------------------------

    assigned_reads_file = (
        input_dir
        / "assigned_reads.csv"
    )


    common_output_dir = (
        input_dir
        / (
            f"common_genes_tp_least"
            f"{args.min_count}_assigned_read_normalized"
        )
    )


    common_output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    common_counts_file = (
        common_output_dir
        / "common_genes_total_counts_summary.csv"
    )


    # --------------------------------------------------------
    # Expected barcodes
    # --------------------------------------------------------

    expected_barcodes = [
        f"bc{i:02d}"
        for i in range(1, 9)
    ]


    print("========================================")
    print("HELIOS NAD-Seq")
    print("Prepare variable 3' normalization inputs")
    print("========================================")
    print(f"Input directory: {input_dir}")
    print(
        f"Time points:     "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )
    print(f"Minimum count:   {args.min_count}")
    print("========================================")


    # ========================================================
    # Read each time point
    # ========================================================

    total_counts_by_tp = {}

    assigned_reads = {}


    for tp_number in range(
        args.start_tp,
        args.end_tp + 1
    ):

        tp_name = f"tp{tp_number}"

        input_file = (
            input_dir
            / tp_name
            / args.input_name
        )


        if not input_file.is_file():

            raise FileNotFoundError(
                f"Missing count matrix for {tp_name}:\n"
                f"{input_file}"
            )


        print()
        print("----------------------------------------")
        print(f"[{tp_name}] Reading:")
        print(f"    {input_file}")


        df = pd.read_csv(
            input_file
        )


        # ----------------------------------------------------
        # Validate Geneid
        # ----------------------------------------------------

        if "Geneid" not in df.columns:

            raise ValueError(
                f"[{tp_name}] 'Geneid' column is missing."
            )


        # ----------------------------------------------------
        # Validate barcodes
        # ----------------------------------------------------

        missing_barcodes = [
            bc
            for bc in expected_barcodes
            if bc not in df.columns
        ]


        if missing_barcodes:

            raise ValueError(
                f"[{tp_name}] Missing barcodes: "
                f"{', '.join(missing_barcodes)}"
            )


        # ----------------------------------------------------
        # Normalize gene IDs
        # ----------------------------------------------------

        df["Geneid"] = (
            df["Geneid"]
            .astype(str)
            .str.strip()
        )


        # ----------------------------------------------------
        # Convert barcode counts to numeric
        # ----------------------------------------------------

        counts = df[
            expected_barcodes
        ].apply(
            pd.to_numeric,
            errors="raise"
        )


        # ----------------------------------------------------
        # Total count per gene for this time point
        #
        # Sum bc01-bc08.
        # ----------------------------------------------------

        gene_totals = counts.sum(
            axis=1
        )


        tp_df = pd.DataFrame(
            {
                "Geneid": df["Geneid"],
                tp_name: gene_totals
            }
        )


        # In case duplicate Geneid entries exist,
        # combine them safely.
        tp_df = (
            tp_df
            .groupby(
                "Geneid",
                as_index=False
            )[tp_name]
            .sum()
        )


        total_counts_by_tp[
            tp_name
        ] = tp_df


        # ----------------------------------------------------
        # Assigned read depth
        #
        # Total number of counts assigned to annotated
        # variable-3' features across all genes and barcodes.
        # ----------------------------------------------------

        assigned_reads[
            tp_name
        ] = int(
            gene_totals.sum()
        )


        print(
            f"[{tp_name}] Genes:          "
            f"{len(tp_df)}"
        )

        print(
            f"[{tp_name}] Assigned reads: "
            f"{assigned_reads[tp_name]}"
        )


    # ========================================================
    # Generate assigned_reads.csv
    # ========================================================

    tp_names = [
        f"tp{i}"
        for i in range(
            args.start_tp,
            args.end_tp + 1
        )
    ]


    assigned_reads_df = pd.DataFrame(
        [
            {
                tp: assigned_reads[tp]
                for tp in tp_names
            }
        ]
    )


    assigned_reads_df.to_csv(
        assigned_reads_file,
        index=False
    )


    print()
    print("========================================")
    print("Assigned-read table generated:")
    print(assigned_reads_file)


    # ========================================================
    # Merge total gene counts across all time points
    # ========================================================

    merged_df = None


    for tp_name in tp_names:

        tp_df = total_counts_by_tp[
            tp_name
        ]


        if merged_df is None:

            merged_df = tp_df.copy()

        else:

            merged_df = pd.merge(
                merged_df,
                tp_df,
                on="Geneid",
                how="inner"
            )


    # --------------------------------------------------------
    # At this point "inner" means that only genes represented
    # at every time point remain.
    # --------------------------------------------------------

    n_common_before_filter = len(
        merged_df
    )


    # --------------------------------------------------------
    # Require >= min-count at EVERY time point
    # --------------------------------------------------------

    keep = (
        merged_df[
            tp_names
        ]
        >= args.min_count
    ).all(
        axis=1
    )


    common_df = merged_df[
        keep
    ].copy()


    common_df = common_df.sort_values(
        "Geneid"
    ).reset_index(
        drop=True
    )


    # --------------------------------------------------------
    # Save common-gene total counts
    # --------------------------------------------------------

    common_df.to_csv(
        common_counts_file,
        index=False
    )


    # ========================================================
    # Summary
    # ========================================================

    print()
    print("========================================")
    print("Preparation completed.")
    print("----------------------------------------")

    print(
        f"Genes present at every tp: "
        f"{n_common_before_filter}"
    )

    print(
        f"Genes >= {args.min_count} reads "
        f"at every tp: {len(common_df)}"
    )

    print()

    print("Assigned reads:")
    print(
        assigned_reads_file
    )

    print()

    print("Common gene counts:")
    print(
        common_counts_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
