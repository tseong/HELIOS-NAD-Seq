#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# SI table preparation
#
# Generate a summary table for NAD-capped genes across tp1-tp16.
#
# Only genes significantly enriched in at least 3 time points
# are retained in the final output.
#
# Final output columns:
#
#   gene_name
#   gene_category
#   NumTimePoints
#   TimePoints
#   Best_adjusted_p_value
#   log2FoldChange_max
#   BaseMean_max
#
# Usage:
#
# python scripts/SI_table_figure_prep/common_nad_genes_with_stats.py \
#     <WORKFLOW_DIR>
#
# Example:
#
# python scripts/SI_table_figure_prep/common_nad_genes_with_stats.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow
#
# Optional:
#   --min-timepoints 3
#   --start-tp 1
#   --end-tp 16
#   --output PATH
# ============================================================


# ------------------------------------------------------------
# Parse GTF attributes
# ------------------------------------------------------------

ATTRIBUTE_RE = re.compile(
    r'(\S+)\s+"([^"]+)"'
)


def parse_gtf_attributes(attribute_string):
    """
    Parse a GTF attribute string into a dictionary.
    """

    return {
        key: value
        for key, value
        in ATTRIBUTE_RE.findall(
            attribute_string
        )
    }


# ------------------------------------------------------------
# Load gene annotation
# ------------------------------------------------------------

def load_gene_annotation(gtf_file):
    """
    Build:
        gene_id -> {
            gene_name,
            gene_category
        }
    """

    annotations = {}


    with gtf_file.open("r") as handle:

        for line in handle:

            if line.startswith("#"):
                continue


            cols = line.rstrip("\n").split("\t")


            if len(cols) < 9:
                continue


            # Only gene-level annotations
            if cols[2] != "gene":
                continue


            attrs = parse_gtf_attributes(
                cols[8]
            )


            gene_id = attrs.get(
                "gene_id"
            )


            if not gene_id:
                continue


            # ------------------------------------------------
            # Gene name
            # ------------------------------------------------

            gene_name = (
                attrs.get("gene")
                or attrs.get("gene_name")
                or attrs.get("locus_tag")
                or gene_id
            )


            # ------------------------------------------------
            # Gene category / biotype
            # ------------------------------------------------

            gene_category = (
                attrs.get("gene_biotype")
                or attrs.get("gene_type")
                or attrs.get("transcript_biotype")
                or attrs.get("transcript_type")
                or "unknown"
            )


            annotations[
                str(gene_id).strip()
            ] = {
                "gene_name": gene_name,
                "gene_category": gene_category,
            }


    return annotations


# ------------------------------------------------------------
# Find column from possible aliases
# ------------------------------------------------------------

def find_column(
    df,
    candidates,
    required=True
):

    lower_map = {
        str(col).lower(): col
        for col in df.columns
    }


    for candidate in candidates:

        if candidate.lower() in lower_map:

            return lower_map[
                candidate.lower()
            ]


    if required:

        raise ValueError(
            "Could not find any of these columns: "
            + ", ".join(candidates)
            + "\nAvailable columns: "
            + ", ".join(
                map(str, df.columns)
            )
        )


    return None


# ------------------------------------------------------------
# Natural timepoint sorting
# ------------------------------------------------------------

def tp_sort_key(tp):

    match = re.search(
        r"(\d+)",
        str(tp)
    )

    if match:

        return int(
            match.group(1)
        )


    return 999999


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare a supplementary-information summary table "
            "for HELIOS NAD genes detected across multiple "
            "time points."
        )
    )


    parser.add_argument(
        "workflow_dir",
        help="Root tp_workflow directory."
    )


    parser.add_argument(
        "--min-timepoints",
        type=int,
        default=3,
        help=(
            "Minimum number of significant time points required "
            "for a gene to be included. Default: 3."
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
        "--output",
        default=None,
        help=(
            "Output CSV. Default: "
            "<workflow_dir>/SI_table_prep/"
            "common_nad_genes_across_timepoints_with_stats.csv"
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Validate arguments
    # --------------------------------------------------------

    if args.min_timepoints < 1:

        raise ValueError(
            "--min-timepoints must be >= 1."
        )


    if args.start_tp > args.end_tp:

        raise ValueError(
            "--start-tp must be <= --end-tp."
        )


    # --------------------------------------------------------
    # Resolve paths
    # --------------------------------------------------------

    workflow_dir = Path(
        args.workflow_dir
    ).expanduser().resolve()


    if not workflow_dir.is_dir():

        raise FileNotFoundError(
            f"Workflow directory does not exist:\n"
            f"{workflow_dir}"
        )


    intergenic_dir = (
        workflow_dir
        / "count_tables"
        / "intergenic"
    )


    common_file = (
        intergenic_dir
        / "common_nad_genes_across_timepoints.csv"
    )


    gtf_file = (
        workflow_dir
        / "gtf"
        / "ncbi_dataset"
        / "data"
        / "GCF_000005845.2"
        / "genomic.gtf"
    )


    if args.output:

        output_file = Path(
            args.output
        ).expanduser().resolve()

    else:

        output_file = (
            workflow_dir
            / "SI_table_prep"
            / "common_nad_genes_across_timepoints_with_stats.csv"
        )


    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Validate files
    # --------------------------------------------------------

    if not common_file.is_file():

        raise FileNotFoundError(
            "Common NAD gene file not found:\n"
            f"{common_file}"
        )


    if not gtf_file.is_file():

        raise FileNotFoundError(
            "GCF annotation GTF not found:\n"
            f"{gtf_file}"
        )


    # --------------------------------------------------------
    # Print configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq")
    print("SI table preparation")
    print("========================================")

    print(
        f"Common NAD genes:\n"
        f"  {common_file}"
    )

    print(
        f"GTF annotation:\n"
        f"  {gtf_file}"
    )

    print(
        f"Time points: "
        f"tp{args.start_tp}-tp{args.end_tp}"
    )

    print(
        f"Minimum required time points: "
        f"{args.min_timepoints}"
    )

    print(
        f"Output:\n"
        f"  {output_file}"
    )

    print("========================================")


    # ========================================================
    # Load common NAD genes
    # ========================================================

    common_df = pd.read_csv(
        common_file
    )


    geneid_col = find_column(
        common_df,
        [
            "Geneid",
            "gene_id",
            "GeneID"
        ]
    )


    common_df[
        geneid_col
    ] = (
        common_df[
            geneid_col
        ]
        .astype(str)
        .str.strip()
    )


    common_genes = list(
        dict.fromkeys(
            common_df[
                geneid_col
            ]
            .dropna()
            .tolist()
        )
    )


    print(
        f"Common NAD genes loaded: "
        f"{len(common_genes)}"
    )


    # ========================================================
    # Load GTF annotation
    # ========================================================

    annotation = load_gene_annotation(
        gtf_file
    )


    print(
        f"Gene annotations loaded: "
        f"{len(annotation)}"
    )


    # ========================================================
    # Read significant results across all time points
    # ========================================================

    stats_by_gene = {
        gene: []
        for gene in common_genes
    }


    common_gene_set = set(
        common_genes
    )


    for tp_num in range(
        args.start_tp,
        args.end_tp + 1
    ):

        tp = f"tp{tp_num}"


        significant_file = (
            intergenic_dir
            / tp
            / "significant_intergenic_Astart.csv"
        )


        if not significant_file.is_file():

            print(
                f"WARNING: Missing {tp} significant file:"
            )

            print(
                f"  {significant_file}"
            )

            continue


        print(
            f"Reading {tp}: "
            f"{significant_file}"
        )


        df = pd.read_csv(
            significant_file
        )


        # ----------------------------------------------------
        # Find relevant columns
        # ----------------------------------------------------

        gene_col = find_column(
            df,
            [
                "Geneid",
                "gene_id",
                "GeneID",
                "Unnamed: 0"
            ]
        )


        padj_col = find_column(
            df,
            [
                "padj",
                "adjusted_p_value",
                "adjusted p-value"
            ]
        )


        log2fc_col = find_column(
            df,
            [
                "log2FoldChange",
                "log2FC",
                "log2foldchange"
            ]
        )


        basemean_col = find_column(
            df,
            [
                "baseMean",
                "basemean"
            ]
        )


        # ----------------------------------------------------
        # Standardize data
        # ----------------------------------------------------

        df[
            gene_col
        ] = (
            df[
                gene_col
            ]
            .astype(str)
            .str.strip()
        )


        df[
            padj_col
        ] = pd.to_numeric(
            df[
                padj_col
            ],
            errors="coerce"
        )


        df[
            log2fc_col
        ] = pd.to_numeric(
            df[
                log2fc_col
            ],
            errors="coerce"
        )


        df[
            basemean_col
        ] = pd.to_numeric(
            df[
                basemean_col
            ],
            errors="coerce"
        )


        # ----------------------------------------------------
        # Keep only common NAD genes
        # ----------------------------------------------------

        df = df[
            df[
                gene_col
            ].isin(
                common_gene_set
            )
        ].copy()


        # ----------------------------------------------------
        # Store time-point-specific statistics
        # ----------------------------------------------------

        for _, row in df.iterrows():

            gene_id = str(
                row[
                    gene_col
                ]
            ).strip()


            stats_by_gene[
                gene_id
            ].append(
                {
                    "timepoint": tp,
                    "padj": row[
                        padj_col
                    ],
                    "log2FoldChange": row[
                        log2fc_col
                    ],
                    "baseMean": row[
                        basemean_col
                    ],
                }
            )


    # ========================================================
    # Build summary
    # ========================================================

    output_rows = []


    for gene_id in common_genes:

        records = stats_by_gene.get(
            gene_id,
            []
        )


        records = sorted(
            records,
            key=lambda x:
            tp_sort_key(
                x[
                    "timepoint"
                ]
            )
        )


        timepoints = [
            record[
                "timepoint"
            ]
            for record in records
        ]


        num_timepoints = len(
            timepoints
        )


        # ----------------------------------------------------
        # Skip genes present in fewer than required TPs
        # ----------------------------------------------------

        if num_timepoints < args.min_timepoints:

            continue


        # ----------------------------------------------------
        # Numeric statistics
        # ----------------------------------------------------

        padj_values = [
            record[
                "padj"
            ]
            for record in records
            if pd.notna(
                record[
                    "padj"
                ]
            )
        ]


        log2fc_values = [
            record[
                "log2FoldChange"
            ]
            for record in records
            if pd.notna(
                record[
                    "log2FoldChange"
                ]
            )
        ]


        basemean_values = [
            record[
                "baseMean"
            ]
            for record in records
            if pd.notna(
                record[
                    "baseMean"
                ]
            )
        ]


        # ----------------------------------------------------
        # Summary statistics
        # ----------------------------------------------------

        best_padj = (
            min(
                padj_values
            )
            if padj_values
            else np.nan
        )


        max_log2fc = (
            max(
                log2fc_values
            )
            if log2fc_values
            else np.nan
        )


        max_basemean = (
            max(
                basemean_values
            )
            if basemean_values
            else np.nan
        )


        # ----------------------------------------------------
        # Annotation
        # ----------------------------------------------------

        meta = annotation.get(
            gene_id,
            {}
        )


        gene_name = meta.get(
            "gene_name",
            gene_id
        )


        gene_category = meta.get(
            "gene_category",
            "unknown"
        )


        # ----------------------------------------------------
        # Final output row
        #
        # Geneid is deliberately NOT included.
        # ----------------------------------------------------

        output_rows.append(
            {
                "gene_name": gene_name,
                "gene_category": gene_category,
                "NumTimePoints": num_timepoints,
                "TimePoints": "; ".join(
                    timepoints
                ),
                "Best_adjusted_p_value": best_padj,
                "log2FoldChange_max": max_log2fc,
                "BaseMean_max": max_basemean,
            }
        )


    # ========================================================
    # Final DataFrame
    # ========================================================

    output_df = pd.DataFrame(
        output_rows,
        columns=[
            "gene_name",
            "gene_category",
            "NumTimePoints",
            "TimePoints",
            "Best_adjusted_p_value",
            "log2FoldChange_max",
            "BaseMean_max",
        ]
    )


    # --------------------------------------------------------
    # Sort:
    # most frequently detected first,
    # then strongest adjusted p-value
    # --------------------------------------------------------

    output_df = output_df.sort_values(
        by=[
            "NumTimePoints",
            "Best_adjusted_p_value"
        ],
        ascending=[
            False,
            True
        ],
        na_position="last"
    ).reset_index(
        drop=True
    )


    # ========================================================
    # Summary
    # ========================================================

    print()
    print("========================================")
    print("Summary")
    print("----------------------------------------")

    print(
        f"Genes with >= "
        f"{args.min_timepoints} time points: "
        f"{len(output_df)}"
    )


    if len(output_df) > 0:

        print(
            f"Maximum NumTimePoints: "
            f"{output_df['NumTimePoints'].max()}"
        )

        print(
            f"Minimum NumTimePoints: "
            f"{output_df['NumTimePoints'].min()}"
        )


    # ========================================================
    # Write output
    # ========================================================

    output_df.to_csv(
        output_file,
        index=False
    )


    print()
    print(
        "Saved SI table:"
    )

    print(
        output_file
    )

    print("========================================")


if __name__ == "__main__":
    main()
