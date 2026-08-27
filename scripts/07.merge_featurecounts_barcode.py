#!/usr/bin/env python3

import pandas as pd
import glob
import os
import re
import argparse


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 07: Merge featureCounts tables
#
# Expected input structure from Step 06:
#
# featurecounts/
# ├── intergenic/
# │   ├── bc01_tp1_eColi_paired.Astart.table
# │   ├── bc01_tp1_eColi_R1_unpaired.Astart.table
# │   └── ...
# │
# └── variable_3prime/
#     ├── bc01_tp1_eColi_paired.Astart.table
#     ├── bc01_tp1_eColi_R1_unpaired.Astart.table
#     └── ...
#
#
# Output structure:
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
#   python scripts/07.merge_counts.py \
#       <FEATURECOUNTS_DIR> \
#       [OUTPUT_DIR]
#
#
# Example:
#
#   python scripts/07.merge_counts.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/featurecounts
#
# ============================================================


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def extract_barcode(filename):
    """
    Extract barcode identifier from filename.

    Example:
        bc01_tp1_eColi_paired.Astart.table
        -> bc01
    """

    match = re.search(
        r'(bc\d+)',
        os.path.basename(filename)
    )

    return match.group(1) if match else None


def extract_timepoint(filename):
    """
    Extract time point from filename.

    Example:
        bc01_tp1_eColi_paired.Astart.table
        -> tp1
    """

    match = re.search(
        r'(tp\d+)',
        os.path.basename(filename)
    )

    return match.group(1) if match else None


def timepoint_sort_key(tp):
    """
    Sort tp1, tp2, ..., tp16 numerically.
    """

    match = re.search(
        r'\d+',
        tp
    )

    return int(match.group()) if match else 0


def barcode_sort_key(barcode):
    """
    Sort bc01, bc02, ..., numerically.
    """

    match = re.search(
        r'\d+',
        barcode
    )

    return int(match.group()) if match else 0


# ------------------------------------------------------------
# Read one featureCounts table
# ------------------------------------------------------------

def read_featurecounts_table(file_path):
    """
    Read Geneid and count column from a featureCounts table.

    featureCounts format:

        Geneid
        Chr
        Start
        End
        Strand
        Length
        <sample count>

    Returns:
        DataFrame with columns:
            Geneid
            count
    """

    df = pd.read_csv(
        file_path,
        sep="\t",
        comment="#",
        header=0
    )


    if df.shape[1] < 7:

        raise ValueError(
            "Unexpected featureCounts format:\n"
            f"{file_path}\n"
            f"Columns detected: {df.shape[1]}"
        )


    # First column = Geneid
    # Last/count column = seventh column for one input SAM

    df = df.iloc[:, [0, 6]].copy()

    df.columns = [
        "Geneid",
        "count"
    ]


    df["count"] = pd.to_numeric(
        df["count"],
        errors="raise"
    )


    return df


# ------------------------------------------------------------
# Merge all tables belonging to one barcode
# ------------------------------------------------------------

def merge_barcode_tables(files, barcode):
    """
    Sum paired + unpaired featureCounts tables belonging to
    a single barcode.
    """

    barcode_df = None


    for file_path in files:

        df = read_featurecounts_table(
            file_path
        )


        if barcode_df is None:

            barcode_df = df

        else:

            barcode_df = pd.merge(
                barcode_df,
                df,
                on="Geneid",
                how="outer",
                suffixes=(
                    "",
                    "_new"
                )
            ).fillna(0)


            barcode_df["count"] = (
                barcode_df["count"]
                + barcode_df["count_new"]
            )


            barcode_df.drop(
                columns=["count_new"],
                inplace=True
            )


    barcode_df.rename(
        columns={
            "count": barcode
        },
        inplace=True
    )


    return barcode_df


# ------------------------------------------------------------
# Process one annotation type
# ------------------------------------------------------------

def process_annotation(
    annotation_name,
    input_dir,
    output_dir
):
    """
    Merge featureCounts tables for one annotation type:
        intergenic
        variable_3prime
    """

    print()
    print("========================================")
    print(f"Annotation: {annotation_name}")
    print("========================================")
    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")


    # --------------------------------------------------------
    # Find all featureCounts tables
    # --------------------------------------------------------

    all_table_files = glob.glob(
        os.path.join(
            input_dir,
            "*.table"
        )
    )


    if not all_table_files:

        print(
            f"WARNING: No .table files found in {input_dir}"
        )

        return 0


    print(
        f"featureCounts tables found: "
        f"{len(all_table_files)}"
    )


    # --------------------------------------------------------
    # Detect time points
    # --------------------------------------------------------

    timepoints = sorted(
        {
            extract_timepoint(file_path)
            for file_path in all_table_files
            if extract_timepoint(file_path) is not None
        },
        key=timepoint_sort_key
    )


    if not timepoints:

        print(
            "WARNING: No time-point identifiers were detected."
        )

        return 0


    print(
        "Time points detected: "
        + ", ".join(timepoints)
    )


    processed_timepoints = 0


    # ========================================================
    # Process each time point
    # ========================================================

    for tp in timepoints:

        print()
        print("----------------------------------------")
        print(f"Processing {annotation_name}: {tp}")


        # ----------------------------------------------------
        # Find all tables for this time point
        # ----------------------------------------------------

        tp_files = glob.glob(
            os.path.join(
                input_dir,
                f"*{tp}*.table"
            )
        )


        # ----------------------------------------------------
        # Paired tables
        #
        # "_unpaired" contains the word "paired", so explicitly
        # exclude unpaired files.
        # ----------------------------------------------------

        paired_files = [
            file_path
            for file_path in glob.glob(
                os.path.join(
                    input_dir,
                    f"*{tp}*_paired*table"
                )
            )
            if "_unpaired" not in os.path.basename(
                file_path
            )
        ]


        # ----------------------------------------------------
        # Unpaired tables
        # ----------------------------------------------------

        unpaired_files = glob.glob(
            os.path.join(
                input_dir,
                f"*{tp}*_unpaired*table"
            )
        )


        print(
            f"  Total tables:    {len(tp_files)}"
        )

        print(
            f"  Paired tables:   {len(paired_files)}"
        )

        print(
            f"  Unpaired tables: {len(unpaired_files)}"
        )


        # ----------------------------------------------------
        # Require paired tables.
        #
        # Unpaired files may legitimately contain zero counts,
        # but their featureCounts tables should still normally
        # exist.
        # ----------------------------------------------------

        if not paired_files:

            print(
                f"WARNING: No paired tables found for {tp}. "
                "Skipping."
            )

            continue


        # ----------------------------------------------------
        # Group by barcode
        # ----------------------------------------------------

        paired_map = {}
        unpaired_map = {}


        for file_path in paired_files:

            barcode = extract_barcode(
                file_path
            )

            if barcode:

                paired_map.setdefault(
                    barcode,
                    []
                ).append(
                    file_path
                )


        for file_path in unpaired_files:

            barcode = extract_barcode(
                file_path
            )

            if barcode:

                unpaired_map.setdefault(
                    barcode,
                    []
                ).append(
                    file_path
                )


        # ----------------------------------------------------
        # Use every barcode found in paired reads.
        #
        # If unpaired counts exist, they are added.
        # ----------------------------------------------------

        barcodes = sorted(
            paired_map.keys(),
            key=barcode_sort_key
        )


        if not barcodes:

            print(
                f"WARNING: No barcode identifiers found "
                f"for {tp}."
            )

            continue


        print(
            "  Barcodes: "
            + ", ".join(barcodes)
        )


        # ----------------------------------------------------
        # Create timepoint output directory
        # ----------------------------------------------------

        tp_output_dir = os.path.join(
            output_dir,
            tp
        )


        os.makedirs(
            tp_output_dir,
            exist_ok=True
        )


        output_filename = os.path.join(
            tp_output_dir,
            "merged_by_barcode_Astart_readCount.csv"
        )


        # ----------------------------------------------------
        # Skip already completed time points
        # ----------------------------------------------------

        if (
            os.path.isfile(output_filename)
            and os.path.getsize(output_filename) > 0
        ):

            print(
                "  Merged output already exists:"
            )

            print(
                f"  {output_filename}"
            )

            print(
                "  Skipping."
            )

            processed_timepoints += 1

            continue


        # ----------------------------------------------------
        # Merge counts
        # ----------------------------------------------------

        merged_df = None


        for barcode in barcodes:

            # Paired table(s)
            all_files = list(
                paired_map.get(
                    barcode,
                    []
                )
            )


            # Add R1/R2 unpaired tables if present
            all_files.extend(
                unpaired_map.get(
                    barcode,
                    []
                )
            )


            print(
                f"  {barcode}: merging "
                f"{len(all_files)} tables"
            )


            barcode_df = merge_barcode_tables(
                all_files,
                barcode
            )


            # ------------------------------------------------
            # Add barcode to merged matrix
            # ------------------------------------------------

            if merged_df is None:

                merged_df = barcode_df

            else:

                merged_df = pd.merge(
                    merged_df,
                    barcode_df,
                    on="Geneid",
                    how="outer"
                ).fillna(0)


        # ----------------------------------------------------
        # Arrange barcode columns
        # ----------------------------------------------------

        barcode_columns = sorted(
            [
                column
                for column in merged_df.columns
                if column.startswith("bc")
            ],
            key=barcode_sort_key
        )


        merged_df = merged_df[
            ["Geneid"]
            + barcode_columns
        ]


        # ----------------------------------------------------
        # Convert counts to integers
        # ----------------------------------------------------

        for barcode in barcode_columns:

            merged_df[barcode] = (
                merged_df[barcode]
                .fillna(0)
                .astype(int)
            )


        # ----------------------------------------------------
        # Write output
        # ----------------------------------------------------

        merged_df.to_csv(
            output_filename,
            index=False
        )


        print(
            "  Merged table written:"
        )

        print(
            f"  {output_filename}"
        )


        print(
            f"  Features: {len(merged_df)}"
        )

        print(
            f"  Barcodes: {len(barcode_columns)}"
        )


        processed_timepoints += 1


    return processed_timepoints


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 07: merge paired and unpaired "
            "featureCounts tables by barcode and time point "
            "for both intergenic and variable 3' annotations."
        )
    )


    parser.add_argument(
        "input_dir",
        help=(
            "Step 06 featureCounts directory containing "
            "'intergenic' and 'variable_3prime' subdirectories."
        )
    )


    parser.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help=(
            "Output directory for merged count matrices. "
            "Default: <parent_of_input_dir>/count_tables"
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Input directory
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


    # --------------------------------------------------------
    # Expected Step 06 directories
    # --------------------------------------------------------

    intergenic_input = os.path.join(
        input_dir,
        "intergenic"
    )


    variable3_input = os.path.join(
        input_dir,
        "variable_3prime"
    )


    if not os.path.isdir(
        intergenic_input
    ):

        raise FileNotFoundError(
            "Intergenic featureCounts directory does not exist:\n"
            f"{intergenic_input}"
        )


    if not os.path.isdir(
        variable3_input
    ):

        raise FileNotFoundError(
            "Variable 3' featureCounts directory does not exist:\n"
            f"{variable3_input}"
        )


    # --------------------------------------------------------
    # Output directory
    # --------------------------------------------------------

    if args.output_dir:

        output_dir = os.path.abspath(
            args.output_dir
        )

    else:

        output_dir = os.path.join(
            os.path.dirname(
                input_dir
            ),
            "count_tables"
        )


    intergenic_output = os.path.join(
        output_dir,
        "intergenic"
    )


    variable3_output = os.path.join(
        output_dir,
        "variable_3prime"
    )


    os.makedirs(
        intergenic_output,
        exist_ok=True
    )


    os.makedirs(
        variable3_output,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Print configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq: Step 07 - Merge counts")
    print("========================================")

    print(
        f"Step 06 directory: {input_dir}"
    )

    print(
        f"Output directory:  {output_dir}"
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
    # Process intergenic counts
    # ========================================================

    n_intergenic = process_annotation(
        annotation_name="intergenic",
        input_dir=intergenic_input,
        output_dir=intergenic_output
    )


    # ========================================================
    # Process variable 3' counts
    # ========================================================

    n_variable3 = process_annotation(
        annotation_name="variable_3prime",
        input_dir=variable3_input,
        output_dir=variable3_output
    )


    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print()
    print("========================================")
    print("Step 07 completed.")
    print("----------------------------------------")

    print(
        f"Intergenic time points:      "
        f"{n_intergenic}"
    )

    print(
        f"Variable 3' time points:     "
        f"{n_variable3}"
    )

    print(
        f"Output directory:"
    )

    print(
        output_dir
    )

    print("========================================")


if __name__ == "__main__":
    main()
