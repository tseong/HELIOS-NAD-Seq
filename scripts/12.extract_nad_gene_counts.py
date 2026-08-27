#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 12: Extract NAD-gene read counts
#
# This script uses:
#
#   common_nad_genes_across_timepoints.csv
#
# generated from Step 11 and extracts the corresponding genes
# from the merged A-start/intergenic read-count matrix at each
# time point.
#
#
# Expected input structure:
#
# count_tables/
# └── intergenic/
#     ├── common_nad_genes_across_timepoints.csv
#     │
#     ├── tp1/
#     │   └── merged_by_barcode_Astart_readCount.csv
#     │
#     ├── tp2/
#     │   └── merged_by_barcode_Astart_readCount.csv
#     │
#     └── ...
#
#
# Output in each tp folder:
#
#   nad_genes_readCount.csv
#
#
# Usage:
#
#   python scripts/12.extract_nad_gene_counts.py \
#       <INTERGENIC_COUNT_DIR>
#
#
# Example:
#
#   python scripts/12.extract_nad_gene_counts.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables/intergenic
#
# ============================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 12: extract read counts for "
            "the global set of NAD genes identified across "
            "time points."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Intergenic count-table directory containing "
            "tp1, tp2, ... subdirectories and "
            "common_nad_genes_across_timepoints.csv."
        )
    )


    parser.add_argument(
        "--nad-gene-file",
        default=None,
        help=(
            "CSV containing the NAD gene list. "
            "Default: "
            "<input_dir>/common_nad_genes_across_timepoints.csv"
        )
    )


    parser.add_argument(
        "--start-tp",
        type=int,
        default=1,
        help=(
            "First time point to process. "
            "Default: 1."
        )
    )


    parser.add_argument(
        "--end-tp",
        type=int,
        default=16,
        help=(
            "Last time point to process. "
            "Default: 16."
        )
    )


    parser.add_argument(
        "--input-name",
        default="merged_by_barcode_Astart_readCount.csv",
        help=(
            "Input count filename inside each tp directory. "
            "Default: merged_by_barcode_Astart_readCount.csv"
        )
    )


    parser.add_argument(
        "--output-name",
        default="nad_genes_readCount.csv",
        help=(
            "Output filename inside each tp directory. "
            "Default: nad_genes_readCount.csv"
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


    if args.nad_gene_file:

        nad_genes_file = Path(
            args.nad_gene_file
        ).expanduser().resolve()

    else:

        nad_genes_file = (
            input_dir
            / "common_nad_genes_across_timepoints.csv"
        )


    if not nad_genes_file.is_file():

        raise FileNotFoundError(
            "NAD gene list does not exist:\n"
            f"{nad_genes_file}"
        )


    # --------------------------------------------------------
    # Validate arguments
    # --------------------------------------------------------

    if args.start_tp < 1:

        raise ValueError(
            "--start-tp must be >= 1."
        )


    if args.end_tp < args.start_tp:

        raise ValueError(
            "--end-tp must be >= --start-tp."
        )


    # --------------------------------------------------------
    # Load NAD gene list
    # --------------------------------------------------------

    nad_df = pd.read_csv(
        nad_genes_file,
        dtype=str
    )


    if "Geneid" not in nad_df.columns:

        raise ValueError(
            "NAD gene list must contain a 'Geneid' column."
        )


    all_genes = (
        nad_df["Geneid"]
        .dropna()
        .astype(str)
        .str.strip()
    )


    all_genes = sorted(
        set(
            gene
            for gene in all_genes
            if gene != ""
        )
    )


    if not all_genes:

        raise ValueError(
            "No valid Geneid values were found in "
            f"{nad_genes_file}"
        )


    # Use a set for fast membership checking
    all_genes_set = set(
        all_genes
    )


    # --------------------------------------------------------
    # Configuration summary
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 12 - Extract NAD gene counts")
    print("========================================")

    print(
        f"Input directory: {input_dir}"
    )

    print(
        f"NAD gene list:   {nad_genes_file}"
    )

    print(
        f"Time points:     "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        f"Target genes:    {len(all_genes)}"
    )

    print(
        f"Input filename:  {args.input_name}"
    )

    print(
        f"Output filename: {args.output_name}"
    )

    print("========================================")


    processed = 0
    skipped = 0

    total_rows_written = 0


    # ========================================================
    # Process each time point
    # ========================================================

    for tp_number in range(
        args.start_tp,
        args.end_tp + 1
    ):

        tp_name = f"tp{tp_number}"

        tp_dir = (
            input_dir
            / tp_name
        )


        # ----------------------------------------------------
        # Check time-point directory
        # ----------------------------------------------------

        if not tp_dir.is_dir():

            print(
                f"[{tp_name}] Directory not found. Skipping."
            )

            skipped += 1
            continue


        input_file = (
            tp_dir
            / args.input_name
        )


        output_file = (
            tp_dir
            / args.output_name
        )


        # ----------------------------------------------------
        # Check count matrix
        # ----------------------------------------------------

        if not input_file.is_file():

            print(
                f"[{tp_name}] Input count file not found. "
                "Skipping:"
            )

            print(
                f"    {input_file}"
            )

            skipped += 1
            continue


        print()
        print("----------------------------------------")

        print(
            f"[{tp_name}] Processing:"
        )

        print(
            f"    {input_file}"
        )


        # ----------------------------------------------------
        # Read count matrix
        # ----------------------------------------------------

        df = pd.read_csv(
            input_file,
            dtype=str
        )


        if "Geneid" not in df.columns:

            print(
                f"[{tp_name}] 'Geneid' column not found. "
                "Skipping."
            )

            skipped += 1
            continue


        # ----------------------------------------------------
        # Normalize Geneid formatting
        # ----------------------------------------------------

        df["Geneid"] = (
            df["Geneid"]
            .astype(str)
            .str.strip()
        )


        # ----------------------------------------------------
        # Extract NAD genes
        # ----------------------------------------------------

        filtered = df[
            df["Geneid"].isin(
                all_genes_set
            )
        ].copy()


        # ----------------------------------------------------
        # Sort consistently by Geneid
        # ----------------------------------------------------

        filtered = filtered.sort_values(
            by="Geneid"
        ).reset_index(
            drop=True
        )


        # ----------------------------------------------------
        # Save output
        # ----------------------------------------------------

        filtered.to_csv(
            output_file,
            index=False
        )


        n_found = len(
            filtered
        )


        n_missing = (
            len(all_genes)
            - n_found
        )


        print(
            f"[{tp_name}] NAD genes found: "
            f"{n_found}/{len(all_genes)}"
        )

        print(
            f"[{tp_name}] NAD genes missing: "
            f"{n_missing}"
        )

        print(
            "[{0}] Output:".format(
                tp_name
            )
        )

        print(
            f"    {output_file}"
        )


        total_rows_written += n_found

        processed += 1


    # ========================================================
    # Final summary
    # ========================================================

    print()
    print("========================================")
    print("Step 12 completed.")
    print("----------------------------------------")

    print(
        f"Processed time points: "
        f"{processed}"
    )

    print(
        f"Skipped time points:   "
        f"{skipped}"
    )

    print(
        f"Target NAD genes:      "
        f"{len(all_genes)}"
    )

    print(
        f"Total rows written:    "
        f"{total_rows_written}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
