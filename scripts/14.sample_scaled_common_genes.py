#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 14: Sample and scale common variable-3' genes
#
# This step:
#
#   1. Loads the Step 13 assigned-read-normalized common-gene
#      count matrix.
#
#   2. Calculates the SUM across tp1-tp16 for each gene.
#
#   3. Sorts genes by SUM in descending order.
#
#   4. Divides genes into three approximately equal abundance
#      bins.
#
#   5. Randomly selects 3 genes from each bin.
#
#   6. Scales bins 2 and 3 so that their median SUM values
#      match the median SUM of bin 1.
#
#   7. Repeats this sampling 10 times.
#
#   8. Generates one summary CSV containing summed tp1-tp16
#      values for each sampling replicate.
#
#
# Expected input:
#
# count_tables/
# └── variable_3prime/
#     └── common_genes_tp_least10_assigned_read_normalized/
#         └── common_genes_total_counts_summary_assigned_reads_normalized.csv
#
#
# Output:
#
# common_genes_tp_least10_assigned_read_normalized/
# └── tp_log2FC_common_genes_assigned_reads_normalized/
#     ├── sampled_1.csv
#     ├── sampled_2.csv
#     ├── ...
#     ├── sampled_10.csv
#     └── sampling_summary.csv
#
#
# Usage:
#
#   python scripts/14.sample_scaled_common_genes.py \
#       <VARIABLE_3PRIME_DIR>
#
#
# Example:
#
#   python scripts/14.sample_scaled_common_genes.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables/variable_3prime
#
# ============================================================


def scale_values(
    df_subset,
    category_median,
    reference_median,
    tp_columns
):
    """
    Scale tp1-tp16 values so that the median SUM of a bin
    matches the median SUM of the reference bin.
    """

    df_subset = df_subset.copy()

    if len(df_subset) == 0:
        return df_subset

    if category_median > 0:

        scaling_factor = (
            reference_median
            / category_median
        )

        df_subset.loc[
            :,
            tp_columns
        ] = (
            df_subset[
                tp_columns
            ]
            * scaling_factor
        )

        df_subset.loc[
            :,
            "SUM"
        ] = (
            df_subset["SUM"]
            * scaling_factor
        )

    return df_subset


def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 14: divide normalized common "
            "variable-3' genes into abundance bins, randomly "
            "sample genes from each bin, scale abundance bins, "
            "and generate replicate summary tables."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Variable 3' count-table directory containing the "
            "Step 13 normalization subdirectory."
        )
    )


    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help=(
            "Minimum-count threshold used in earlier steps. "
            "Used to identify the normalization directory. "
            "Default: 10."
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
        "--genes-per-bin",
        type=int,
        default=3,
        help=(
            "Number of genes sampled from each abundance bin. "
            "Default: 3."
        )
    )


    parser.add_argument(
        "--n-samplings",
        type=int,
        default=10,
        help=(
            "Number of independent sampling replicates. "
            "Default: 10."
        )
    )


    parser.add_argument(
        "--seed-start",
        type=int,
        default=0,
        help=(
            "Starting random seed. Sampling replicate i uses "
            "seed-start + i. Default: 0."
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


    if args.genes_per_bin < 1:

        raise ValueError(
            "--genes-per-bin must be >= 1."
        )


    if args.n_samplings < 1:

        raise ValueError(
            "--n-samplings must be >= 1."
        )


    tp_columns = [
        f"tp{i}"
        for i in range(
            args.start_tp,
            args.end_tp + 1
        )
    ]


    # --------------------------------------------------------
    # Step 13 input
    # --------------------------------------------------------

    normalization_dir = (
        input_dir
        / (
            f"common_genes_tp_least"
            f"{args.min_count}_assigned_read_normalized"
        )
    )


    input_file = (
        normalization_dir
        / "common_genes_total_counts_summary_assigned_reads_normalized.csv"
    )


    if not input_file.is_file():

        raise FileNotFoundError(
            "Step 13 normalized count file does not exist:\n"
            f"{input_file}"
        )


    # --------------------------------------------------------
    # Output directory
    # --------------------------------------------------------

    output_dir = (
        normalization_dir
        / "tp_log2FC_common_genes_assigned_reads_normalized"
    )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    summary_output_file = (
        output_dir
        / "sampling_summary.csv"
    )


    # --------------------------------------------------------
    # Configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 14 - Sample scaled common genes")
    print("========================================")

    print(
        f"Input file:       {input_file}"
    )

    print(
        f"Output directory: {output_dir}"
    )

    print(
        f"Time points:      "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        f"Genes per bin:    {args.genes_per_bin}"
    )

    print(
        f"Sampling runs:    {args.n_samplings}"
    )

    print("========================================")


    # ========================================================
    # Load normalized counts
    # ========================================================

    df = pd.read_csv(
        input_file
    )


    # --------------------------------------------------------
    # Validate time-point columns
    # --------------------------------------------------------

    missing_tp = [
        tp
        for tp in tp_columns
        if tp not in df.columns
    ]


    if missing_tp:

        raise ValueError(
            "Input file is missing time-point columns: "
            + ", ".join(missing_tp)
        )


    # --------------------------------------------------------
    # Convert time-point data to numeric
    # --------------------------------------------------------

    df[
        tp_columns
    ] = df[
        tp_columns
    ].apply(
        pd.to_numeric,
        errors="raise"
    )


    # ========================================================
    # Calculate total abundance
    # ========================================================

    df["SUM"] = df[
        tp_columns
    ].sum(
        axis=1
    )


    # Sort highest-abundance genes first
    df = df.sort_values(
        by="SUM",
        ascending=False
    ).reset_index(
        drop=True
    )


    # ========================================================
    # Divide into three bins
    # ========================================================

    total_rows = len(
        df
    )


    if total_rows == 0:

        raise ValueError(
            "No genes were found in the normalized input file."
        )


    bin_size, remainder = divmod(
        total_rows,
        3
    )


    num_bin1 = bin_size
    num_bin2 = bin_size
    num_bin3 = bin_size


    # Preserve original distribution strategy
    if remainder == 1:

        num_bin2 += 1

    elif remainder == 2:

        num_bin2 += 1
        num_bin3 += 1


    bin1 = df.iloc[
        :num_bin1
    ].copy()


    bin2 = df.iloc[
        num_bin1:
        num_bin1 + num_bin2
    ].copy()


    bin3 = df.iloc[
        num_bin1 + num_bin2:
    ].copy()


    # --------------------------------------------------------
    # Median SUM values
    # --------------------------------------------------------

    median_bin1 = bin1[
        "SUM"
    ].median()


    median_bin2 = bin2[
        "SUM"
    ].median()


    median_bin3 = bin3[
        "SUM"
    ].median()


    print()
    print("Abundance bins:")
    print(
        f"  Bin 1: {len(bin1)} genes, "
        f"median SUM = {median_bin1}"
    )

    print(
        f"  Bin 2: {len(bin2)} genes, "
        f"median SUM = {median_bin2}"
    )

    print(
        f"  Bin 3: {len(bin3)} genes, "
        f"median SUM = {median_bin3}"
    )


    if len(bin1) == 0:

        raise ValueError(
            "Bin 1 contains no genes."
        )


    if median_bin1 <= 0:

        raise ValueError(
            "Median SUM of Bin 1 must be > 0."
        )


    # ========================================================
    # Generate sampling replicates
    # ========================================================

    summary_data = []


    for sampling_index in range(
        args.n_samplings
    ):

        random_seed = (
            args.seed_start
            + sampling_index
        )


        # ----------------------------------------------------
        # Sample Bin 1
        # ----------------------------------------------------

        if len(bin1) > 0:

            sample_bin1 = bin1.sample(
                n=min(
                    len(bin1),
                    args.genes_per_bin
                ),
                random_state=random_seed,
                replace=False
            ).copy()

        else:

            sample_bin1 = pd.DataFrame(
                columns=df.columns
            )


        # ----------------------------------------------------
        # Sample and scale Bin 2
        # ----------------------------------------------------

        if len(bin2) > 0:

            sample_bin2 = bin2.sample(
                n=min(
                    len(bin2),
                    args.genes_per_bin
                ),
                random_state=random_seed,
                replace=False
            ).copy()


            sample_bin2 = scale_values(
                sample_bin2,
                category_median=median_bin2,
                reference_median=median_bin1,
                tp_columns=tp_columns
            )

        else:

            sample_bin2 = pd.DataFrame(
                columns=df.columns
            )


        # ----------------------------------------------------
        # Sample and scale Bin 3
        # ----------------------------------------------------

        if len(bin3) > 0:

            sample_bin3 = bin3.sample(
                n=min(
                    len(bin3),
                    args.genes_per_bin
                ),
                random_state=random_seed,
                replace=False
            ).copy()


            sample_bin3 = scale_values(
                sample_bin3,
                category_median=median_bin3,
                reference_median=median_bin1,
                tp_columns=tp_columns
            )

        else:

            sample_bin3 = pd.DataFrame(
                columns=df.columns
            )


        # ----------------------------------------------------
        # Combine sampled genes
        # ----------------------------------------------------

        sampled_df = pd.concat(
            [
                sample_bin1,
                sample_bin2,
                sample_bin3
            ],
            ignore_index=True
        )


        # ----------------------------------------------------
        # Save replicate
        # ----------------------------------------------------

        sampling_number = (
            sampling_index + 1
        )


        output_file = (
            output_dir
            / f"sampled_{sampling_number}.csv"
        )


        sampled_df.to_csv(
            output_file,
            index=False
        )


        print(
            f"Sampling {sampling_number}: "
            f"{len(sampled_df)} genes -> "
            f"{output_file}"
        )


        # ----------------------------------------------------
        # Summary totals
        # ----------------------------------------------------

        tp_sums = (
            sampled_df[
                tp_columns
            ]
            .sum(
                axis=0
            )
        )


        total_sum = (
            sampled_df[
                "SUM"
            ]
            .sum()
        )


        summary_row = {
            "Sampling": (
                f"Sampling_{sampling_number}"
            )
        }


        for tp in tp_columns:

            summary_row[
                tp
            ] = tp_sums[
                tp
            ]


        summary_row[
            "SUM"
        ] = total_sum


        summary_data.append(
            summary_row
        )


    # ========================================================
    # Generate summary
    # ========================================================

    summary_df = pd.DataFrame(
        summary_data,
        columns=[
            "Sampling",
            *tp_columns,
            "SUM"
        ]
    )


    summary_df.to_csv(
        summary_output_file,
        index=False
    )


    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print()
    print("========================================")
    print("Step 14 completed.")
    print("----------------------------------------")

    print(
        f"Total input genes: "
        f"{total_rows}"
    )

    print(
        f"Sampling replicates: "
        f"{args.n_samplings}"
    )

    print(
        f"Genes requested per bin: "
        f"{args.genes_per_bin}"
    )

    print()

    print(
        "Summary output:"
    )

    print(
        summary_output_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
