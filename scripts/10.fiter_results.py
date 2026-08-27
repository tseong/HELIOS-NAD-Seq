#!/usr/bin/env python3

import argparse
import os
import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 10: Filter PyDESeq2 results
#
# Expected input structure:
#
# count_tables/
# ├── intergenic/
# │   ├── tp1/
# │   │   └── *_results.csv
# │   └── ...
# │
# └── variable_3prime/
#     ├── tp1/
#     │   └── *_results.csv
#     └── ...
#
#
# Filtering:
#
# intergenic:
#   padj < 0.05
#   log2FoldChange >= 1
#
# variable_3prime:
#   -1 <= log2FoldChange <= 1
#   baseMean >= 100
#
#
# Usage:
#
#   python scripts/10.filter_results.py \
#       <COUNT_TABLES_DIR>
#
# Example:
#
#   python scripts/10.filter_results.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables
#
# ============================================================


def find_results_file(tp_dir, annotation_name):
    """
    Find the single Step 09 *_results.csv file in a time-point
    directory.
    """

    result_files = sorted(
        [
            filename
            for filename in os.listdir(tp_dir)
            if filename.endswith("_results.csv")
        ]
    )

    if not result_files:

        return None

    if len(result_files) > 1:

        raise RuntimeError(
            f"[{annotation_name}] More than one *_results.csv "
            f"file found in {tp_dir}: {result_files}"
        )

    return os.path.join(
        tp_dir,
        result_files[0]
    )


def filter_intergenic(
    df,
    padj_threshold,
    log2fc_threshold
):
    """
    Filter A-start/intergenic DESeq2 results.

    Criteria:
        padj < threshold
        log2FoldChange >= threshold
    """

    required_columns = {
        "log2FoldChange",
        "padj"
    }

    missing = (
        required_columns
        - set(df.columns)
    )

    if missing:

        raise ValueError(
            "Intergenic results are missing required columns: "
            + ", ".join(sorted(missing))
        )


    valid_df = df.dropna(
        subset=[
            "padj",
            "log2FoldChange"
        ]
    )


    filtered_df = valid_df[
        (valid_df["padj"] < padj_threshold)
        &
        (
            valid_df["log2FoldChange"]
            >= log2fc_threshold
        )
    ].copy()


    filtered_df = filtered_df.sort_values(
        by=[
            "padj",
            "log2FoldChange"
        ],
        ascending=[
            True,
            False
        ]
    )


    return filtered_df, len(valid_df)


def filter_variable3(
    df,
    min_log2fc,
    max_log2fc,
    min_basemean
):
    """
    Filter variable 3' DESeq2 results.

    Criteria:
        min_log2fc <= log2FoldChange <= max_log2fc
        baseMean >= min_basemean

    padj is intentionally not used for this filter.
    """

    required_columns = {
        "log2FoldChange",
        "baseMean"
    }

    missing = (
        required_columns
        - set(df.columns)
    )

    if missing:

        raise ValueError(
            "Variable 3' results are missing required columns: "
            + ", ".join(sorted(missing))
        )


    valid_df = df.dropna(
        subset=[
            "log2FoldChange",
            "baseMean"
        ]
    )


    filtered_df = valid_df[
        (
            valid_df["log2FoldChange"]
            >= min_log2fc
        )
        &
        (
            valid_df["log2FoldChange"]
            <= max_log2fc
        )
        &
        (
            valid_df["baseMean"]
            >= min_basemean
        )
    ].copy()


    filtered_df = filtered_df.sort_values(
        by=[
            "baseMean",
            "log2FoldChange"
        ],
        ascending=[
            False,
            False
        ]
    )


    return filtered_df, len(valid_df)


def process_annotation(
    annotation_name,
    annotation_dir,
    start_tp,
    end_tp,
    args
):

    print()
    print("========================================")
    print(f"Annotation: {annotation_name}")
    print("========================================")

    processed = 0
    skipped = 0
    total_filtered = 0


    for tp in range(
        start_tp,
        end_tp + 1
    ):

        tp_name = f"tp{tp}"

        tp_dir = os.path.join(
            annotation_dir,
            tp_name
        )


        if not os.path.isdir(tp_dir):

            print(
                f"[{annotation_name} / {tp_name}] "
                "Directory not found. Skipping."
            )

            skipped += 1
            continue


        results_path = find_results_file(
            tp_dir,
            annotation_name
        )


        if results_path is None:

            print(
                f"[{annotation_name} / {tp_name}] "
                "No *_results.csv file found. Skipping."
            )

            skipped += 1
            continue


        print()
        print("----------------------------------------")
        print(
            f"[{annotation_name} / {tp_name}] Processing:"
        )
        print(
            f"    {results_path}"
        )


        df = pd.read_csv(
            results_path,
            index_col=0
        )


        # ====================================================
        # Intergenic / A-start branch
        # ====================================================

        if annotation_name == "intergenic":

            filtered_df, n_valid = filter_intergenic(
                df,
                padj_threshold=args.padj,
                log2fc_threshold=args.log2fc
            )


            output_name = (
                "significant_intergenic_Astart.csv"
            )


        # ====================================================
        # Variable 3' branch
        # ====================================================

        elif annotation_name == "variable_3prime":

            filtered_df, n_valid = filter_variable3(
                df,
                min_log2fc=args.variable_min_log2fc,
                max_log2fc=args.variable_max_log2fc,
                min_basemean=args.variable_min_basemean
            )


            output_name = (
                "variable3prime_broad_filtered.csv"
            )


        else:

            raise ValueError(
                f"Unknown annotation type: "
                f"{annotation_name}"
            )


        output_path = os.path.join(
            tp_dir,
            output_name
        )


        filtered_df.to_csv(
            output_path,
            index=True
        )


        n_total = len(df)
        n_filtered = len(filtered_df)


        print(
            f"Total tested features:    {n_total}"
        )

        print(
            f"Features with valid data: {n_valid}"
        )

        print(
            f"Features retained:        {n_filtered}"
        )

        print(
            "Saved:"
        )

        print(
            f"    {output_path}"
        )


        total_filtered += n_filtered
        processed += 1


    print()
    print("----------------------------------------")
    print(
        f"{annotation_name} summary"
    )
    print("----------------------------------------")

    print(
        f"Processed time points: {processed}"
    )

    print(
        f"Skipped time points:   {skipped}"
    )

    print(
        f"Total retained hits:   {total_filtered}"
    )


    return (
        processed,
        skipped,
        total_filtered
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 10: filter PyDESeq2 results "
            "for both intergenic/A-start and variable 3' "
            "analysis branches."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Step 09 count_tables directory containing "
            "'intergenic' and 'variable_3prime' subdirectories."
        )
    )


    parser.add_argument(
        "--start-tp",
        type=int,
        default=1,
        help="First time point to process. Default: 1."
    )


    parser.add_argument(
        "--end-tp",
        type=int,
        default=16,
        help="Last time point to process. Default: 16."
    )


    # --------------------------------------------------------
    # Intergenic filter parameters
    # --------------------------------------------------------

    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help=(
            "Adjusted p-value threshold for intergenic/A-start "
            "results. Default: 0.05."
        )
    )


    parser.add_argument(
        "--log2fc",
        type=float,
        default=1.0,
        help=(
            "Minimum log2FoldChange for intergenic/A-start "
            "results. Default: 1."
        )
    )


    # --------------------------------------------------------
    # Variable 3' filter parameters
    # --------------------------------------------------------

    parser.add_argument(
        "--variable-min-log2fc",
        type=float,
        default=-1.0,
        help=(
            "Minimum log2FoldChange for variable 3' filter. "
            "Default: -1."
        )
    )


    parser.add_argument(
        "--variable-max-log2fc",
        type=float,
        default=1.0,
        help=(
            "Maximum log2FoldChange for variable 3' filter. "
            "Default: 1."
        )
    )


    parser.add_argument(
        "--variable-min-basemean",
        type=float,
        default=100.0,
        help=(
            "Minimum baseMean for variable 3' filter. "
            "Default: 100."
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Validate input
    # --------------------------------------------------------

    input_dir = os.path.abspath(
        args.input_dir
    )


    if not os.path.isdir(
        input_dir
    ):

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


    if not 0 < args.padj <= 1:

        raise ValueError(
            "--padj must be > 0 and <= 1."
        )


    if (
        args.variable_min_log2fc
        > args.variable_max_log2fc
    ):

        raise ValueError(
            "--variable-min-log2fc cannot be greater than "
            "--variable-max-log2fc."
        )


    if args.variable_min_basemean < 0:

        raise ValueError(
            "--variable-min-basemean must be >= 0."
        )


    # --------------------------------------------------------
    # Annotation directories
    # --------------------------------------------------------

    intergenic_dir = os.path.join(
        input_dir,
        "intergenic"
    )


    variable3_dir = os.path.join(
        input_dir,
        "variable_3prime"
    )


    if not os.path.isdir(
        intergenic_dir
    ):

        raise FileNotFoundError(
            "Intergenic directory does not exist:\n"
            f"{intergenic_dir}"
        )


    if not os.path.isdir(
        variable3_dir
    ):

        raise FileNotFoundError(
            "Variable 3' directory does not exist:\n"
            f"{variable3_dir}"
        )


    # --------------------------------------------------------
    # Configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq: Step 10 - Filter results")
    print("========================================")

    print(
        f"Input directory: {input_dir}"
    )

    print(
        f"Time points:     "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print()

    print("Intergenic / A-start filter:")

    print(
        f"  padj < {args.padj}"
    )

    print(
        f"  log2FoldChange >= {args.log2fc}"
    )

    print()

    print("Variable 3' filter:")

    print(
        f"  {args.variable_min_log2fc} "
        "<= log2FoldChange <= "
        f"{args.variable_max_log2fc}"
    )

    print(
        f"  baseMean >= "
        f"{args.variable_min_basemean}"
    )

    print("========================================")


    # ========================================================
    # Intergenic
    # ========================================================

    (
        intergenic_processed,
        intergenic_skipped,
        intergenic_total
    ) = process_annotation(
        annotation_name="intergenic",
        annotation_dir=intergenic_dir,
        start_tp=args.start_tp,
        end_tp=args.end_tp,
        args=args
    )


    # ========================================================
    # Variable 3'
    # ========================================================

    (
        variable_processed,
        variable_skipped,
        variable_total
    ) = process_annotation(
        annotation_name="variable_3prime",
        annotation_dir=variable3_dir,
        start_tp=args.start_tp,
        end_tp=args.end_tp,
        args=args
    )


    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print()
    print("========================================")
    print("Step 10 completed.")
    print("----------------------------------------")

    print(
        f"Intergenic processed:      "
        f"{intergenic_processed}"
    )

    print(
        f"Intergenic skipped:        "
        f"{intergenic_skipped}"
    )

    print(
        f"Intergenic retained hits:  "
        f"{intergenic_total}"
    )

    print()

    print(
        f"Variable 3' processed:     "
        f"{variable_processed}"
    )

    print(
        f"Variable 3' skipped:       "
        f"{variable_skipped}"
    )

    print(
        f"Variable 3' retained hits: "
        f"{variable_total}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
