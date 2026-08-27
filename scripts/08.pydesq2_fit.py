#!/usr/bin/env python3

import argparse
import os
import pickle as pkl

import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 08: Fit PyDESeq2 models
#
# Expected input structure from Step 07:
#
# count_tables/
# ├── intergenic/
# │   ├── tp1/
# │   │   └── merged_by_barcode_Astart_readCount.csv
# │   ├── tp2/
# │   └── ...
# │
# └── variable_3prime/
#     ├── tp1/
#     │   └── merged_by_barcode_Astart_readCount.csv
#     ├── tp2/
#     └── ...
#
#
# Usage:
#
#   python scripts/08.pydeseq2_dds.py \
#       <COUNT_TABLES_DIR>
#
#
# Example:
#
#   python scripts/08.pydeseq2_dds.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables
#
# ============================================================


def process_annotation(
    annotation_name,
    annotation_dir,
    start_tp,
    end_tp,
    threads,
    min_total_count
):
    """
    Fit PyDESeq2 models independently for each time point
    within one annotation branch.
    """

    print()
    print("========================================")
    print(f"Annotation: {annotation_name}")
    print("========================================")
    print(f"Input directory:     {annotation_dir}")
    print(f"Time points:         tp{start_tp}-tp{end_tp}")
    print(f"CPUs:                {threads}")
    print(f"Minimum total count: {min_total_count}")
    print("========================================")


    inference = DefaultInference(
        n_cpus=threads
    )


    processed = 0
    skipped = 0


    for tp in range(
        start_tp,
        end_tp + 1
    ):

        tp_name = f"tp{tp}"

        tp_dir = os.path.join(
            annotation_dir,
            tp_name
        )


        # ----------------------------------------------------
        # Input/output paths
        # ----------------------------------------------------

        input_csv = os.path.join(
            tp_dir,
            "merged_by_barcode_Astart_readCount.csv"
        )


        output_pkl = os.path.join(
            tp_dir,
            "merged_by_barcode_Astart_readCount_dds.pkl"
        )


        # ----------------------------------------------------
        # Skip missing time-point directory
        # ----------------------------------------------------

        if not os.path.isdir(tp_dir):

            print(
                f"[{annotation_name} / {tp_name}] "
                "Directory not found. Skipping."
            )

            skipped += 1
            continue


        # ----------------------------------------------------
        # Skip missing input CSV
        # ----------------------------------------------------

        if not os.path.isfile(input_csv):

            print(
                f"[{annotation_name} / {tp_name}] "
                "Input file not found. Skipping:"
            )

            print(
                f"    {input_csv}"
            )

            skipped += 1
            continue


        # ----------------------------------------------------
        # Skip already fitted DDS
        # ----------------------------------------------------

        if (
            os.path.isfile(output_pkl)
            and os.path.getsize(output_pkl) > 0
        ):

            print(
                f"[{annotation_name} / {tp_name}] "
                "DDS already exists. Skipping:"
            )

            print(
                f"    {output_pkl}"
            )

            processed += 1
            continue


        print()
        print("----------------------------------------")
        print(
            f"[{annotation_name} / {tp_name}] Processing"
        )
        print(
            f"Input: {input_csv}"
        )


        # ----------------------------------------------------
        # Read count matrix
        # ----------------------------------------------------

        df = pd.read_csv(
            input_csv,
            index_col=0
        )


        # ----------------------------------------------------
        # Expected HELIOS barcode structure
        #
        # bc01-bc04 = positive samples
        # bc05-bc08 = negative controls
        # ----------------------------------------------------

        expected_barcodes = [
            f"bc{i:02d}"
            for i in range(1, 9)
        ]


        missing_barcodes = [
            barcode
            for barcode in expected_barcodes
            if barcode not in df.columns
        ]


        if missing_barcodes:

            raise ValueError(
                f"[{annotation_name} / {tp_name}] "
                "Missing expected barcodes: "
                f"{', '.join(missing_barcodes)}"
            )


        # Keep only bc01-bc08 in correct order
        df = df[
            expected_barcodes
        ]


        # ----------------------------------------------------
        # Convert counts to numeric
        # ----------------------------------------------------

        counts = df.apply(
            pd.to_numeric,
            errors="coerce"
        )


        if counts.isna().any().any():

            raise ValueError(
                f"[{annotation_name} / {tp_name}] "
                "Non-numeric or missing values detected "
                "in count matrix."
            )


        # featureCounts produces integer read counts
        counts = counts.astype(int)


        # ----------------------------------------------------
        # Check for negative counts
        # ----------------------------------------------------

        if (counts < 0).any().any():

            raise ValueError(
                f"[{annotation_name} / {tp_name}] "
                "Negative count values detected."
            )


        # ----------------------------------------------------
        # Filter low-count features
        # ----------------------------------------------------

        n_before = counts.shape[0]


        keep = (
            counts.sum(axis=1)
            >= min_total_count
        )


        counts = counts.loc[
            keep
        ]


        n_after = counts.shape[0]


        print(
            f"Features before filtering: {n_before}"
        )

        print(
            f"Features retained:         {n_after}"
        )


        if n_after == 0:

            print(
                f"[{annotation_name} / {tp_name}] "
                "No features remain after filtering. "
                "Skipping."
            )

            skipped += 1
            continue


        # ----------------------------------------------------
        # Metadata
        #
        # bc01-bc04 = HELIOS positive samples
        # bc05-bc08 = negative controls
        #
        # Keep original contrast naming:
        #
        # Treated vs Control
        # ----------------------------------------------------

        metadata = pd.DataFrame(
            {
                "conditions": (
                    ["Treated"] * 4
                    + ["Control"] * 4
                )
            },
            index=expected_barcodes
        )


        # ----------------------------------------------------
        # PyDESeq2 expects:
        #
        # samples x features
        # ----------------------------------------------------

        counts = counts.T


        # ----------------------------------------------------
        # Confirm sample ordering
        # ----------------------------------------------------

        if not counts.index.equals(
            metadata.index
        ):

            raise RuntimeError(
                f"[{annotation_name} / {tp_name}] "
                "Count matrix and metadata sample order "
                "do not match."
            )


        # ----------------------------------------------------
        # Fit DESeq2 model
        # ----------------------------------------------------

        print(
            f"[{annotation_name} / {tp_name}] "
            "Creating DeseqDataSet..."
        )


        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors="conditions",
            refit_cooks=False,
            inference=inference
        )


        print(
            f"[{annotation_name} / {tp_name}] "
            "Fitting size factors..."
        )


        dds.fit_size_factors(
            fit_type="ratio"
        )


        print(
            f"[{annotation_name} / {tp_name}] "
            "Fitting gene-wise dispersions..."
        )

        dds.fit_genewise_dispersions()


        print(
            f"[{annotation_name} / {tp_name}] "
            "Fitting dispersion trend..."
        )

        dds.fit_dispersion_trend()


        print(
            f"[{annotation_name} / {tp_name}] "
            "Fitting dispersion prior..."
        )

        dds.fit_dispersion_prior()


        print(
            f"[{annotation_name} / {tp_name}] "
            "Fitting MAP dispersions..."
        )

        dds.fit_MAP_dispersions()


        print(
            f"[{annotation_name} / {tp_name}] "
            "Fitting log2 fold changes..."
        )

        dds.fit_LFC()


        # ----------------------------------------------------
        # Cook's distance handling
        #
        # refit_cooks=False above, so this normally does not
        # run. Retained for compatibility with the previous
        # workflow.
        # ----------------------------------------------------

        if dds.refit_cooks:

            dds.varm["replaced"] = np.zeros_like(
                dds.var_names,
                dtype=bool
            )

            dds.refit()


        # ----------------------------------------------------
        # Save fitted DDS
        # ----------------------------------------------------

        tmp_pkl = (
            output_pkl
            + ".tmp"
        )


        with open(
            tmp_pkl,
            "wb"
        ) as f:

            pkl.dump(
                dds,
                f
            )


        # Only move into final location after successful save
        os.replace(
            tmp_pkl,
            output_pkl
        )


        print(
            f"[{annotation_name} / {tp_name}] "
            "Saved DDS:"
        )

        print(
            f"    {output_pkl}"
        )


        processed += 1


    print()
    print("----------------------------------------")
    print(
        f"{annotation_name} summary"
    )
    print("----------------------------------------")
    print(
        f"Processed/existing time points: "
        f"{processed}"
    )
    print(
        f"Skipped time points:            "
        f"{skipped}"
    )


    return processed, skipped


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 08: fit PyDESeq2 models "
            "independently for each time point for both "
            "intergenic and variable 3' count matrices."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Step 07 count_tables directory containing "
            "'intergenic' and 'variable_3prime' subdirectories."
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
        "--threads",
        type=int,
        default=8,
        help=(
            "Number of CPUs used by PyDESeq2. "
            "Default: 8."
        )
    )


    parser.add_argument(
        "--min-total-count",
        type=int,
        default=10,
        help=(
            "Minimum total count across all 8 barcode samples "
            "required to retain a feature. Default: 10."
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Validate general arguments
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


    if args.threads < 1:

        raise ValueError(
            "--threads must be >= 1."
        )


    if args.min_total_count < 0:

        raise ValueError(
            "--min-total-count must be >= 0."
        )


    # --------------------------------------------------------
    # Expected Step 07 annotation directories
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
            "Intergenic count-table directory does not exist:\n"
            f"{intergenic_dir}"
        )


    if not os.path.isdir(
        variable3_dir
    ):

        raise FileNotFoundError(
            "Variable 3' count-table directory does not exist:\n"
            f"{variable3_dir}"
        )


    # --------------------------------------------------------
    # Print configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq: Step 08 - PyDESeq2 DDS")
    print("========================================")

    print(
        f"Count-table directory: {input_dir}"
    )

    print(
        f"Time points:           "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        f"CPUs:                  {args.threads}"
    )

    print(
        f"Minimum total count:   "
        f"{args.min_total_count}"
    )

    print(
        "Annotations:"
    )

    print(
        "  intergenic"
    )

    print(
        "  variable_3prime"
    )

    print("========================================")


    # ========================================================
    # Intergenic
    # ========================================================

    intergenic_processed, intergenic_skipped = (
        process_annotation(
            annotation_name="intergenic",
            annotation_dir=intergenic_dir,
            start_tp=args.start_tp,
            end_tp=args.end_tp,
            threads=args.threads,
            min_total_count=args.min_total_count
        )
    )


    # ========================================================
    # Variable 3'
    # ========================================================

    variable3_processed, variable3_skipped = (
        process_annotation(
            annotation_name="variable_3prime",
            annotation_dir=variable3_dir,
            start_tp=args.start_tp,
            end_tp=args.end_tp,
            threads=args.threads,
            min_total_count=args.min_total_count
        )
    )


    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print()
    print("========================================")
    print("Step 08 completed.")
    print("----------------------------------------")

    print(
        f"Intergenic processed/existing: "
        f"{intergenic_processed}"
    )

    print(
        f"Intergenic skipped:            "
        f"{intergenic_skipped}"
    )

    print(
        f"Variable 3' processed/existing:"
        f" {variable3_processed}"
    )

    print(
        f"Variable 3' skipped:           "
        f"{variable3_skipped}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
