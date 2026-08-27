#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 15_3: Normalize NAD-gene counts
#
# Performs:
#   1. Total read-depth normalization
#   2. Common-gene sampling normalization
#
# Usage:
#
#   python scripts/15_3.normalize_nad_gene_counts.py \
#       <MAIN_NAD_TABLE> \
#       <WORKFLOW_DIR>
#
# Example:
#
#   cd /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow
#
#   python scripts/15_3.normalize_nad_gene_counts.py \
#       nad_genes_across_tp_with_readCount_helios.tsv \
#       .
#
# Output columns only:
#
#   gene_name
#   timepoint
#   gene_biotype
#   3PAB_rep1
#   3PAB_rep2
#   3PAB_rep3
#   3PAB_rep4
#
# ============================================================


REPS = [
    "3PAB_rep1",
    "3PAB_rep2",
    "3PAB_rep3",
    "3PAB_rep4",
]


OUTPUT_COLUMNS = [
    "gene_name",
    "timepoint",
    "gene_biotype",
    *REPS,
]


EXPECTED_MAIN_COLS = set(
    OUTPUT_COLUMNS
)


# ============================================================
# Load main NAD table
# ============================================================

def load_main_table(path):

    df = pd.read_csv(
        path,
        sep="\t"
    )


    if not EXPECTED_MAIN_COLS.issubset(
        df.columns
    ):

        df = pd.read_csv(
            path,
            sep=r"\s+",
            engine="python"
        )


    missing = (
        EXPECTED_MAIN_COLS
        - set(df.columns)
    )


    if missing:

        raise ValueError(
            "Main NAD table is missing columns: "
            + ", ".join(
                sorted(missing)
            )
        )


    # Standardize timepoint labels
    df["timepoint"] = (
        df["timepoint"]
        .astype(str)
        .str.strip()
    )


    df["timepoint"] = df[
        "timepoint"
    ].apply(
        lambda x:
        x if x.startswith("tp")
        else f"tp{x}"
    )


    # Force counts to float
    for rep in REPS:

        df[rep] = pd.to_numeric(
            df[rep],
            errors="raise"
        ).astype("float64")


    return df


# ============================================================
# Load one-row tp1-tp16 normalization table
# ============================================================

def load_tp_values(
    path,
    tp_columns
):

    df = pd.read_csv(
        path
    )


    if len(df) == 0:

        raise ValueError(
            f"No rows found in:\n{path}"
        )


    missing = [
        tp
        for tp in tp_columns
        if tp not in df.columns
    ]


    if missing:

        raise ValueError(
            f"{path} is missing timepoints: "
            + ", ".join(missing)
        )


    values = pd.to_numeric(
        df.loc[
            0,
            tp_columns
        ],
        errors="raise"
    ).astype(
        "float64"
    )


    if (
        values <= 0
    ).any():

        invalid = list(
            values[
                values <= 0
            ].index
        )

        raise ValueError(
            "Non-positive normalization values for: "
            + ", ".join(invalid)
        )


    return values


# ============================================================
# Map factor to rows according to time point
# ============================================================

def map_factor(
    df,
    factors,
    factor_name
):

    factor_dict = {
        str(tp): float(value)
        for tp, value
        in factors.items()
    }


    mapped = df[
        "timepoint"
    ].map(
        factor_dict
    )


    if mapped.isna().any():

        bad_tp = sorted(
            set(
                df.loc[
                    mapped.isna(),
                    "timepoint"
                ]
            )
        )

        raise ValueError(
            f"Could not map {factor_name} for: "
            + ", ".join(bad_tp)
        )


    return mapped.astype(
        "float64"
    )


# ============================================================
# Read-depth normalization
# ============================================================

def make_read_depth_normalized(
    raw_df,
    read_depth_factors
):

    result = raw_df.copy()


    # Keep raw values internally for verification
    for rep in REPS:

        result[
            f"{rep}_raw"
        ] = result[
            rep
        ].copy()


    result[
        "read_depth_factor"
    ] = map_factor(
        result,
        read_depth_factors,
        "read-depth factor"
    )


    for rep in REPS:

        result[
            rep
        ] = (
            result[
                f"{rep}_raw"
            ]
            * result[
                "read_depth_factor"
            ]
        )


    return result


# ============================================================
# Read-depth + sampling normalization
# ============================================================

def make_sampling_normalized(
    raw_df,
    read_depth_factors,
    sampling_factors
):

    result = raw_df.copy()


    for rep in REPS:

        result[
            f"{rep}_raw"
        ] = result[
            rep
        ].copy()


    result[
        "read_depth_factor"
    ] = map_factor(
        result,
        read_depth_factors,
        "read-depth factor"
    )


    result[
        "sampling_factor"
    ] = map_factor(
        result,
        sampling_factors,
        "sampling factor"
    )


    result[
        "final_factor"
    ] = (
        result[
            "read_depth_factor"
        ]
        * result[
            "sampling_factor"
        ]
    )


    for rep in REPS:

        result[
            rep
        ] = (
            result[
                f"{rep}_raw"
            ]
            * result[
                "final_factor"
            ]
        )


    return result


# ============================================================
# Verify normalized values
# ============================================================

def verify_output(
    df,
    factor_column,
    label
):

    print()
    print("========================================")
    print(f"VERIFYING {label}")
    print("========================================")


    example = (
        df[
            [
                "timepoint",
                factor_column,
                "3PAB_rep1_raw",
                "3PAB_rep1",
            ]
        ]
        .groupby(
            "timepoint",
            sort=False
        )
        .head(1)
        .copy()
    )


    example["_tp"] = (
        example[
            "timepoint"
        ]
        .str.replace(
            "tp",
            "",
            regex=False
        )
        .astype(int)
    )


    example = (
        example
        .sort_values("_tp")
        .drop(
            columns="_tp"
        )
    )


    print(
        example.to_string(
            index=False
        )
    )


    for rep in REPS:

        expected = (
            df[
                f"{rep}_raw"
            ]
            * df[
                factor_column
            ]
        )


        difference = (
            df[rep]
            - expected
        ).abs()


        max_difference = float(
            difference.max()
        )


        if max_difference > 1e-8:

            raise RuntimeError(
                f"{label}: verification FAILED "
                f"for {rep}. "
                f"Maximum difference = "
                f"{max_difference}"
            )


    print()
    print(
        "Mathematical verification PASSED."
    )


# ============================================================
# Save only final requested columns
# ============================================================

def save_clean_output(
    df,
    output_file
):

    clean_df = df[
        OUTPUT_COLUMNS
    ].copy()


    clean_df.to_csv(
        output_file,
        sep="\t",
        index=False,
        float_format="%.8f"
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 15_3: apply read-depth "
            "and common-gene sampling normalization."
        )
    )


    parser.add_argument(
        "main_file",
        help=(
            "Main NAD table generated by Step 15_2."
        )
    )


    parser.add_argument(
        "workflow_dir",
        help=(
            "Root tp_workflow directory."
        )
    )


    parser.add_argument(
        "--start-tp",
        type=int,
        default=1
    )


    parser.add_argument(
        "--end-tp",
        type=int,
        default=16
    )


    parser.add_argument(
        "--reference-tp",
        default="tp1"
    )


    parser.add_argument(
        "--min-count",
        type=int,
        default=10
    )


    parser.add_argument(
        "--output-dir",
        default=None
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Paths
    # --------------------------------------------------------

    main_file = Path(
        args.main_file
    ).expanduser().resolve()


    workflow_dir = Path(
        args.workflow_dir
    ).expanduser().resolve()


    if not main_file.is_file():

        raise FileNotFoundError(
            f"Main NAD table not found:\n"
            f"{main_file}"
        )


    if not workflow_dir.is_dir():

        raise FileNotFoundError(
            f"Workflow directory not found:\n"
            f"{workflow_dir}"
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
            f"Invalid reference timepoint: "
            f"{args.reference_tp}"
        )


    # --------------------------------------------------------
    # Step 15_1 read depth file
    # --------------------------------------------------------

    read_depth_file = (
        workflow_dir
        / "trimmomatic"
        / "read_depth_by_timepoint_all_reads.csv"
    )


    # --------------------------------------------------------
    # Step 14 sampling summary
    # --------------------------------------------------------

    sampling_file = (
        workflow_dir
        / "count_tables"
        / "variable_3prime"
        / (
            f"common_genes_tp_least"
            f"{args.min_count}_assigned_read_normalized"
        )
        / "tp_log2FC_common_genes_assigned_reads_normalized"
        / "sampling_summary.csv"
    )


    # --------------------------------------------------------
    # Output directory
    # --------------------------------------------------------

    if args.output_dir:

        output_dir = Path(
            args.output_dir
        ).expanduser().resolve()

    else:

        output_dir = (
            workflow_dir
            / "normalized_nad_genes"
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    if not read_depth_file.is_file():

        raise FileNotFoundError(
            f"Read-depth file not found:\n"
            f"{read_depth_file}"
        )


    if not sampling_file.is_file():

        raise FileNotFoundError(
            f"Sampling summary not found:\n"
            f"{sampling_file}"
        )


    # ========================================================
    # Load raw NAD table
    # ========================================================

    raw_df = load_main_table(
        main_file
    )


    print("========================================")
    print("HELIOS NAD-Seq Step 15_3")
    print("========================================")

    print(
        f"Rows in raw NAD table: "
        f"{len(raw_df)}"
    )

    print(
        f"Reference: "
        f"{args.reference_tp}"
    )


    # ========================================================
    # Read-depth factors
    # ========================================================

    read_depths = load_tp_values(
        read_depth_file,
        tp_columns
    )


    reference_depth = float(
        read_depths[
            args.reference_tp
        ]
    )


    read_depth_factors = (
        reference_depth
        / read_depths
    ).astype(
        "float64"
    )


    print()
    print("Read-depth factors:")
    print("----------------------------------------")


    for tp in tp_columns:

        print(
            f"{tp}: "
            f"depth={read_depths[tp]:.0f}, "
            f"factor="
            f"{read_depth_factors[tp]:.8f}"
        )


    # ========================================================
    # Read-depth normalized output
    # ========================================================

    read_depth_df = make_read_depth_normalized(
        raw_df,
        read_depth_factors
    )


    verify_output(
        read_depth_df,
        "read_depth_factor",
        "read-depth normalization"
    )


    read_depth_output = (
        output_dir
        / "read_depth_normalized.tsv"
    )


    save_clean_output(
        read_depth_df,
        read_depth_output
    )


    print()
    print(
        f"Wrote:\n"
        f"{read_depth_output}"
    )


    # ========================================================
    # Load sampling summary
    # ========================================================

    sampling_df = pd.read_csv(
        sampling_file
    )


    required_sampling = {
        "Sampling",
        *tp_columns
    }


    missing = (
        required_sampling
        - set(
            sampling_df.columns
        )
    )


    if missing:

        raise ValueError(
            "Sampling summary is missing: "
            + ", ".join(
                sorted(missing)
            )
        )


    # ========================================================
    # Sampling normalizations
    # ========================================================

    for _, row in sampling_df.iterrows():

        sampling_name = str(
            row[
                "Sampling"
            ]
        ).strip()


        sampling_index = int(
            sampling_name.split(
                "_"
            )[-1]
        )


        sampling_values = pd.to_numeric(
            row[
                tp_columns
            ],
            errors="raise"
        ).astype(
            "float64"
        )


        if (
            sampling_values <= 0
        ).any():

            raise ValueError(
                f"{sampling_name} contains "
                "non-positive values."
            )


        sampling_reference = float(
            sampling_values[
                args.reference_tp
            ]
        )


        sampling_factors = (
            sampling_reference
            / sampling_values
        ).astype(
            "float64"
        )


        print()
        print("========================================")

        print(
            sampling_name
        )

        print("========================================")


        for tp in tp_columns:

            print(
                f"{tp}: "
                f"sampling_factor="
                f"{sampling_factors[tp]:.8f}"
            )


        final_df = make_sampling_normalized(
            raw_df,
            read_depth_factors,
            sampling_factors
        )


        verify_output(
            final_df,
            "final_factor",
            sampling_name
        )


        output_file = (
            output_dir
            / (
                "read_depth_common_genes_normalized_"
                f"{sampling_index}.tsv"
            )
        )


        save_clean_output(
            final_df,
            output_file
        )


        print()
        print(
            f"Wrote:\n"
            f"{output_file}"
        )


    # ========================================================
    # Finish
    # ========================================================

    print()
    print("========================================")
    print("Step 15_3 completed.")
    print("========================================")

    print(
        f"Read-depth normalized:\n"
        f"{read_depth_output}"
    )

    print()

    print(
        f"Sampling-normalized files: "
        f"{len(sampling_df)}"
    )


if __name__ == "__main__":
    main()
