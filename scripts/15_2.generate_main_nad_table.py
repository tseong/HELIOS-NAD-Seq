#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 15_2: Generate main NAD-gene read-count table
#
# This step prepares the main input table required for
# downstream normalization.
#
#
# Input from Step 12:
#
# count_tables/
# └── intergenic/
#     ├── tp1/
#     │   └── nad_genes_readCount.csv
#     ├── tp2/
#     │   └── nad_genes_readCount.csv
#     └── ...
#
#
# GTF:
#
# gtf/
# └── ncbi_dataset/
#     └── data/
#         └── GCF_000005845.2/
#             └── genomic.gtf
#
#
# Output:
#
# nad_genes_across_tp_with_readCount_helios.tsv
#
#
# Output columns:
#
# gene_name
# timepoint
# gene_biotype
# 3PAB_rep1
# 3PAB_rep2
# 3PAB_rep3
# 3PAB_rep4
#
#
# Mapping:
#
# bc01 -> 3PAB_rep1
# bc02 -> 3PAB_rep2
# bc03 -> 3PAB_rep3
# bc04 -> 3PAB_rep4
#
#
# Usage:
#
# python scripts/15_2.generate_main_nad_table.py \
#     <WORKFLOW_DIR>
#
#
# Example:
#
# python scripts/15_2.generate_main_nad_table.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow
#
# ============================================================


# ------------------------------------------------------------
# Parse GTF attributes
# ------------------------------------------------------------

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attrs(attribute_string):
    """
    Parse GTF attributes into a dictionary.
    """

    return {
        key: value
        for key, value in ATTR_RE.findall(
            attribute_string
        )
    }


# ------------------------------------------------------------
# Read gene metadata from GTF
# ------------------------------------------------------------

def load_gene_metadata(gtf_file):
    """
    Build a mapping from gene_id to:

        gene_name
        gene_biotype
    """

    metadata = {}


    with gtf_file.open("r") as handle:

        for line in handle:

            if line.startswith("#"):
                continue


            columns = line.rstrip("\n").split("\t")


            if len(columns) < 9:
                continue


            if columns[2] != "gene":
                continue


            attrs = parse_attrs(
                columns[8]
            )


            gene_id = attrs.get(
                "gene_id"
            )


            if not gene_id:
                continue


            # ------------------------------------------------
            # Gene name
            #
            # Prefer gene_name, then gene, then locus_tag,
            # finally gene_id itself.
            # ------------------------------------------------

            gene_name = (
                attrs.get("gene_name")
                or attrs.get("gene")
                or attrs.get("locus_tag")
                or gene_id
            )


            # ------------------------------------------------
            # Gene biotype
            # ------------------------------------------------

            gene_biotype = (
                attrs.get("gene_biotype")
                or attrs.get("transcript_biotype")
                or attrs.get("gene_type")
                or attrs.get("transcript_type")
                or "unknown"
            )


            metadata[
                gene_id
            ] = {
                "gene_name": gene_name,
                "gene_biotype": gene_biotype
            }


    return metadata


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "HELIOS NAD-Seq Step 15_2: generate the main "
            "NAD-gene read-count table across tp1-tp16 "
            "using Step 12 NAD-gene count tables."
        )
    )


    parser.add_argument(
        "workflow_dir",
        help="Root HELIOS workflow directory."
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
            "Output TSV. Default: "
            "<workflow_dir>/"
            "nad_genes_across_tp_with_readCount_helios.tsv"
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Resolve workflow directory
    # --------------------------------------------------------

    workflow_dir = Path(
        args.workflow_dir
    ).expanduser().resolve()


    if not workflow_dir.is_dir():

        raise FileNotFoundError(
            f"Workflow directory does not exist:\n"
            f"{workflow_dir}"
        )


    if args.start_tp < 1:

        raise ValueError(
            "--start-tp must be >= 1."
        )


    if args.end_tp < args.start_tp:

        raise ValueError(
            "--end-tp must be >= --start-tp."
        )


    # --------------------------------------------------------
    # Input locations
    # --------------------------------------------------------

    intergenic_dir = (
        workflow_dir
        / "count_tables"
        / "intergenic"
    )


    gtf_file = (
        workflow_dir
        / "gtf"
        / "ncbi_dataset"
        / "data"
        / "GCF_000005845.2"
        / "genomic.gtf"
    )


    # --------------------------------------------------------
    # Output
    # --------------------------------------------------------

    if args.output:

        output_file = Path(
            args.output
        ).expanduser().resolve()

    else:

        output_file = (
            workflow_dir
            / "nad_genes_across_tp_with_readCount_helios.tsv"
        )


    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Validate required inputs
    # --------------------------------------------------------

    if not intergenic_dir.is_dir():

        raise FileNotFoundError(
            "Intergenic count-table directory does not exist:\n"
            f"{intergenic_dir}"
        )


    if not gtf_file.is_file():

        raise FileNotFoundError(
            "GCF reference GTF does not exist:\n"
            f"{gtf_file}"
        )


    # --------------------------------------------------------
    # Positive HELIOS barcodes
    # --------------------------------------------------------

    barcode_map = {
        "bc01": "3PAB_rep1",
        "bc02": "3PAB_rep2",
        "bc03": "3PAB_rep3",
        "bc04": "3PAB_rep4",
    }


    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 15_2 - Generate main NAD table")
    print("========================================")

    print(
        f"Workflow directory: {workflow_dir}"
    )

    print(
        f"Intergenic counts:  {intergenic_dir}"
    )

    print(
        f"GCF annotation:     {gtf_file}"
    )

    print(
        f"Time points:        "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        f"Output:             {output_file}"
    )

    print("========================================")


    # ========================================================
    # Load gene metadata
    # ========================================================

    gene_metadata = load_gene_metadata(
        gtf_file
    )


    print(
        f"Genes loaded from GTF: "
        f"{len(gene_metadata)}"
    )


    # ========================================================
    # Process time points
    # ========================================================

    all_rows = []

    processed = 0
    skipped = 0

    metadata_missing = set()


    for tp_number in range(
        args.start_tp,
        args.end_tp + 1
    ):

        tp_name = f"tp{tp_number}"


        input_file = (
            intergenic_dir
            / tp_name
            / "nad_genes_readCount.csv"
        )


        if not input_file.is_file():

            print(
                f"[{tp_name}] "
                "nad_genes_readCount.csv not found. Skipping:"
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
        # Read Step 12 count table
        # ----------------------------------------------------

        df = pd.read_csv(
            input_file
        )


        required_columns = {
            "Geneid",
            *barcode_map.keys()
        }


        missing_columns = (
            required_columns
            - set(df.columns)
        )


        if missing_columns:

            raise ValueError(
                f"[{tp_name}] Missing required columns: "
                + ", ".join(
                    sorted(
                        missing_columns
                    )
                )
            )


        # ----------------------------------------------------
        # Convert barcode counts to numeric
        # ----------------------------------------------------

        for barcode in barcode_map:

            df[
                barcode
            ] = pd.to_numeric(
                df[
                    barcode
                ],
                errors="raise"
            )


        # ----------------------------------------------------
        # Generate output row for each NAD gene
        # ----------------------------------------------------

        for _, row in df.iterrows():

            gene_id = str(
                row[
                    "Geneid"
                ]
            ).strip()


            meta = gene_metadata.get(
                gene_id
            )


            if meta is None:

                metadata_missing.add(
                    gene_id
                )

                gene_name = gene_id
                gene_biotype = "unknown"

            else:

                gene_name = meta[
                    "gene_name"
                ]

                gene_biotype = meta[
                    "gene_biotype"
                ]


            output_row = {
                "gene_name": gene_name,
                "timepoint": tp_name,
                "gene_biotype": gene_biotype,
            }


            for barcode, rep_name in barcode_map.items():

                output_row[
                    rep_name
                ] = row[
                    barcode
                ]


            all_rows.append(
                output_row
            )


        print(
            f"[{tp_name}] NAD genes added: "
            f"{len(df)}"
        )


        processed += 1


    # ========================================================
    # Build final table
    # ========================================================

    if not all_rows:

        raise RuntimeError(
            "No NAD-gene rows were generated."
        )


    output_df = pd.DataFrame(
        all_rows,
        columns=[
            "gene_name",
            "timepoint",
            "gene_biotype",
            "3PAB_rep1",
            "3PAB_rep2",
            "3PAB_rep3",
            "3PAB_rep4",
        ]
    )


    # --------------------------------------------------------
    # Sort consistently
    # --------------------------------------------------------

    output_df[
        "_tp_number"
    ] = (
        output_df[
            "timepoint"
        ]
        .str.replace(
            "tp",
            "",
            regex=False
        )
        .astype(int)
    )


    output_df = output_df.sort_values(
        by=[
            "_tp_number",
            "gene_name"
        ]
    ).drop(
        columns=[
            "_tp_number"
        ]
    ).reset_index(
        drop=True
    )


    # ========================================================
    # Save TSV
    # ========================================================

    output_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )


    # ========================================================
    # Summary
    # ========================================================

    print()
    print("========================================")
    print("Step 15_2 completed.")
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
        f"Total NAD-gene rows:   "
        f"{len(output_df)}"
    )

    print(
        f"Unique gene names:     "
        f"{output_df['gene_name'].nunique()}"
    )

    print(
        f"Genes without GTF metadata: "
        f"{len(metadata_missing)}"
    )


    if metadata_missing:

        print()
        print(
            "WARNING: Gene IDs without matching GTF metadata:"
        )

        for gene_id in sorted(
            metadata_missing
        ):

            print(
                f"  {gene_id}"
            )


    print()
    print(
        "Main NAD table:"
    )

    print(
        output_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
