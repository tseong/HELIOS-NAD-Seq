#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Generate common_nad_genes_across_timepoints.csv
#
# This script collects significantly enriched NAD genes from
# the Step 11 intergenic/A-start results across all time points.
#
# Expected directory structure:
#
# count_tables/
# └── intergenic/
#     ├── tp1/
#     │   └── significant_intergenic_Astart.csv
#     ├── tp2/
#     │   └── significant_intergenic_Astart.csv
#     ├── ...
#     └── tp16/
#         └── significant_intergenic_Astart.csv
#
#
# Output:
#
# common_nad_genes_across_timepoints.csv
#
#
# Usage:
#
#   python scripts/11.generate_common_nad_genes.py \
#       <INTERGENIC_COUNT_DIR>
#
#
# Example:
#
#   python scripts/11.generate_common_nad_genes.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables/intergenic
#
# ============================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate the union of significant HELIOS NAD genes "
            "identified from significant_intergenic_Astart.csv "
            "across multiple time points."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Intergenic count-table directory containing "
            "tp1, tp2, ... subdirectories."
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
        "--output",
        default=None,
        help=(
            "Output CSV file. Default: "
            "<input_dir>/common_nad_genes_across_timepoints.csv"
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Resolve input directory
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


    # --------------------------------------------------------
    # Output path
    # --------------------------------------------------------

    if args.output:

        output_file = Path(
            args.output
        ).expanduser().resolve()

    else:

        output_file = (
            input_dir
            / "common_nad_genes_across_timepoints.csv"
        )


    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq")
    print("Generate common NAD gene list")
    print("========================================")

    print(
        f"Input directory: {input_dir}"
    )

    print(
        f"Time points:     "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        "Input filename:  "
        "significant_intergenic_Astart.csv"
    )

    print(
        f"Output:          {output_file}"
    )

    print("========================================")


    # --------------------------------------------------------
    # Store:
    #
    # Geneid -> set of timepoint numbers
    #
    # Example:
    #
    # b0001 -> {1, 2, 5}
    # --------------------------------------------------------

    gene_timepoints = {}


    processed = 0
    skipped = 0
    total_gene_timepoint_hits = 0


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
        # Check tp directory
        # ----------------------------------------------------

        if not tp_dir.is_dir():

            print(
                f"[{tp_name}] Directory not found. Skipping."
            )

            skipped += 1

            continue


        # ----------------------------------------------------
        # Step 10 intergenic significant file
        # ----------------------------------------------------

        input_file = (
            tp_dir
            / "significant_intergenic_Astart.csv"
        )


        if not input_file.is_file():

            print(
                f"[{tp_name}] File not found. Skipping:"
            )

            print(
                f"    {input_file}"
            )

            skipped += 1

            continue


        print()
        print("----------------------------------------")
        print(
            f"[{tp_name}] Reading:"
        )

        print(
            f"    {input_file}"
        )


        # ----------------------------------------------------
        # Read significant gene table
        # ----------------------------------------------------

        df = pd.read_csv(
            input_file,
            dtype=str
        )


        # ----------------------------------------------------
        # Determine Geneid column
        #
        # Depending on how the results CSV was saved, Geneid
        # may either be explicitly named "Geneid" or stored as
        # the first index column.
        # ----------------------------------------------------

        if "Geneid" in df.columns:

            gene_column = "Geneid"

        elif len(df.columns) > 0:

            gene_column = df.columns[0]

        else:

            print(
                f"[{tp_name}] Empty CSV with no columns. "
                "Skipping."
            )

            skipped += 1

            continue


        # ----------------------------------------------------
        # Extract gene IDs
        # ----------------------------------------------------

        genes = (
            df[gene_column]
            .dropna()
            .astype(str)
            .str.strip()
        )


        genes = genes[
            genes != ""
        ]


        # Remove duplicates within one time point
        genes = set(
            genes
        )


        print(
            f"[{tp_name}] Significant NAD genes: "
            f"{len(genes)}"
        )


        # ----------------------------------------------------
        # Add to global gene collection
        # ----------------------------------------------------

        for gene_id in genes:

            if gene_id not in gene_timepoints:

                gene_timepoints[
                    gene_id
                ] = set()


            gene_timepoints[
                gene_id
            ].add(
                tp_number
            )


        total_gene_timepoint_hits += len(
            genes
        )

        processed += 1


    # ========================================================
    # Build output table
    # ========================================================

    records = []


    for gene_id, timepoints in gene_timepoints.items():

        sorted_timepoints = sorted(
            timepoints
        )


        timepoint_names = [
            f"tp{tp}"
            for tp in sorted_timepoints
        ]


        records.append(
            {
                "Geneid": gene_id,
                "n_timepoints": len(
                    sorted_timepoints
                ),
                "timepoints": ",".join(
                    timepoint_names
                )
            }
        )


    common_df = pd.DataFrame(
        records,
        columns=[
            "Geneid",
            "n_timepoints",
            "timepoints"
        ]
    )


    # --------------------------------------------------------
    # Sort genes:
    #
    # 1. Most frequently detected genes first
    # 2. Geneid alphabetically
    # --------------------------------------------------------

    if not common_df.empty:

        common_df = common_df.sort_values(
            by=[
                "n_timepoints",
                "Geneid"
            ],
            ascending=[
                False,
                True
            ]
        ).reset_index(
            drop=True
        )


    # --------------------------------------------------------
    # Write result
    # --------------------------------------------------------

    common_df.to_csv(
        output_file,
        index=False
    )


    # ========================================================
    # Summary
    # ========================================================

    print()
    print("========================================")
    print("Common NAD gene list generated")
    print("----------------------------------------")

    print(
        f"Processed time points:      "
        f"{processed}"
    )

    print(
        f"Skipped time points:        "
        f"{skipped}"
    )

    print(
        f"Unique NAD genes:           "
        f"{len(common_df)}"
    )

    print(
        f"Gene-timepoint detections:  "
        f"{total_gene_timepoint_hits}"
    )


    # --------------------------------------------------------
    # Frequency distribution
    # --------------------------------------------------------

    if not common_df.empty:

        print()
        print("Detection frequency:")


        frequency = (
            common_df[
                "n_timepoints"
            ]
            .value_counts()
            .sort_index(
                ascending=False
            )
        )


        for n_timepoints, n_genes in frequency.items():

            print(
                f"  {int(n_timepoints):2d} "
                f"time points: "
                f"{n_genes} genes"
            )


    print()
    print("Output:")
    print(
        output_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
