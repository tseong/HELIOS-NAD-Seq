#!/usr/bin/env python3

import argparse
import math
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import hypergeom


# ============================================================
# HELIOS NAD-Seq pipeline
# SI figure/table preparation
#
# KEGG pathway enrichment for:
#
#   1. DTW Cluster 1 genes
#   2. DTW Cluster 2 genes
#   3. tp7 significant NAD genes
#
# Script location:
#
#   scripts/SI_figure_table_prep/
#   17.kegg_enrichment_clusters_tp7.py
#
# Usage:
#
#   python scripts/SI_figure_table_prep/17.kegg_enrichment_clusters_tp7.py \
#       <WORKFLOW_DIR>
#
# Example:
#
#   python scripts/SI_figure_table_prep/17.kegg_enrichment_clusters_tp7.py \
#       /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow
#
# Optional:
#
#   --cluster-dir PATH
#   --timepoint 7
#   --top-n 10
#   --padj 0.05
#   --output-dir PATH
#
# Default output:
#
#   <WORKFLOW_DIR>/SI_figure_table_prep/kegg_enrichment/
#
# ============================================================n 10
#   --padj 0.05
#   --output-dir PATH
#
# ============================================================


# ============================================================
# GTF parsing
# ============================================================

ATTR_RE = re.compile(
    r'(\S+)\s+"([^"]+)"'
)


def parse_gtf_attributes(text):

    return {
        key: value
        for key, value
        in ATTR_RE.findall(text)
    }


# ============================================================
# Gene identifier normalization
# ============================================================

def normalize_gene_id(value):

    value = str(value).strip()


    if ":" in value:

        value = value.split(":")[-1]


    for prefix in [
        "gene-",
        "rna-",
        "cds-",
    ]:

        if value.startswith(prefix):

            value = value[
                len(prefix):
            ]


    return value


# ============================================================
# Load GTF genes
# ============================================================

def load_gtf_genes(gtf_file):
    """
    Returns:
        set of normalized E. coli gene/locus identifiers
    """

    genes = set()


    with gtf_file.open("r") as handle:

        for line in handle:

            if line.startswith("#"):
                continue


            cols = line.rstrip(
                "\n"
            ).split(
                "\t"
            )


            if len(cols) < 9:
                continue


            if cols[2] != "gene":
                continue


            attrs = parse_gtf_attributes(
                cols[8]
            )


            gene_id = attrs.get(
                "gene_id"
            )


            locus_tag = attrs.get(
                "locus_tag"
            )


            if gene_id:

                genes.add(
                    normalize_gene_id(
                        gene_id
                    )
                )


            if locus_tag:

                genes.add(
                    normalize_gene_id(
                        locus_tag
                    )
                )


    return genes


# ============================================================
# KEGG download/cache
# ============================================================

def get_kegg_text(
    url,
    cache_file
):

    if cache_file.is_file():

        print(
            f"Using cached KEGG file: "
            f"{cache_file}"
        )

        return cache_file.read_text()


    print(
        f"Downloading KEGG data:"
    )

    print(
        f"  {url}"
    )


    try:

        with urllib.request.urlopen(
            url,
            timeout=60
        ) as response:

            text = (
                response
                .read()
                .decode(
                    "utf-8"
                )
            )

    except Exception as exc:

        raise RuntimeError(
            "Could not download KEGG data and no cache "
            f"was available:\n{cache_file}\n\n"
            f"{exc}"
        )


    cache_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    cache_file.write_text(
        text
    )


    return text


# ============================================================
# Clean KEGG pathway names
# ============================================================

def clean_pathway_name(name):

    name = str(
        name
    ).strip()


    # Remove:
    #
    # " - Escherichia coli K-12 MG1655"
    #
    name = re.sub(
        r"\s*-\s*Escherichia coli K-12 MG1655\s*$",
        "",
        name
    )


    return name.strip()


# ============================================================
# Load KEGG pathway names
# ============================================================

def load_kegg_pathway_names(
    cache_dir
):

    text = get_kegg_text(
        "https://rest.kegg.jp/list/pathway/eco",
        cache_dir
        / "kegg_ecoli_pathways.txt"
    )


    pathway_names = {}


    for line in text.splitlines():

        if not line.strip():
            continue


        parts = line.split(
            "\t",
            1
        )


        if len(parts) != 2:
            continue


        pathway_id = (
            parts[0]
            .replace(
                "path:",
                ""
            )
            .strip()
        )


        pathway_name = clean_pathway_name(
            parts[1]
        )


        pathway_names[
            pathway_id
        ] = pathway_name


    return pathway_names


# ============================================================
# Load KEGG gene-pathway associations
# ============================================================

def load_kegg_gene_pathways(
    cache_dir,
    pathway_names
):

    text = get_kegg_text(
        "https://rest.kegg.jp/link/pathway/eco",
        cache_dir
        / "kegg_ecoli_gene_pathway_links.txt"
    )


    gene_to_pathways = {}

    pathway_to_genes = {}


    for line in text.splitlines():

        if not line.strip():
            continue


        parts = line.split(
            "\t"
        )


        if len(parts) != 2:
            continue


        gene = normalize_gene_id(
            parts[0]
        )


        pathway_id = (
            parts[1]
            .replace(
                "path:",
                ""
            )
            .strip()
        )


        pathway_name = pathway_names.get(
            pathway_id,
            pathway_id
        )


        gene_to_pathways.setdefault(
            gene,
            set()
        ).add(
            pathway_name
        )


        pathway_to_genes.setdefault(
            pathway_name,
            set()
        ).add(
            gene
        )


    return (
        gene_to_pathways,
        pathway_to_genes
    )


# ============================================================
# Benjamini-Hochberg correction
# ============================================================

def benjamini_hochberg(
    pvalues
):

    pvalues = np.asarray(
        pvalues,
        dtype=float
    )


    n = len(
        pvalues
    )


    if n == 0:

        return np.array(
            []
        )


    order = np.argsort(
        pvalues
    )


    ranked = pvalues[
        order
    ]


    adjusted_ranked = (
        ranked
        * n
        / np.arange(
            1,
            n + 1
        )
    )


    # Enforce monotonicity
    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[
            ::-1
        ]
    )[
        ::-1
    ]


    adjusted_ranked = np.clip(
        adjusted_ranked,
        0,
        1
    )


    adjusted = np.empty(
        n,
        dtype=float
    )


    adjusted[
        order
    ] = adjusted_ranked


    return adjusted


# ============================================================
# Hypergeometric KEGG enrichment
# ============================================================

def run_kegg_enrichment(
    query_genes,
    background_genes,
    pathway_to_genes
):

    query_genes = set(
        query_genes
    )


    background_genes = set(
        background_genes
    )


    # Only genes in background are eligible
    query_genes = (
        query_genes
        & background_genes
    )


    M = len(
        background_genes
    )


    N = len(
        query_genes
    )


    if N == 0:

        raise ValueError(
            "No query genes overlap with the KEGG background."
        )


    results = []


    for pathway, pathway_genes in (
        pathway_to_genes.items()
    ):

        pathway_background = (
            pathway_genes
            & background_genes
        )


        n = len(
            pathway_background
        )


        overlap = (
            query_genes
            & pathway_background
        )


        k = len(
            overlap
        )


        if k == 0:
            continue


        # ----------------------------------------------------
        # Hypergeometric:
        #
        # probability of >= k genes occurring in pathway
        # ----------------------------------------------------

        pvalue = hypergeom.sf(
            k - 1,
            M,
            n,
            N
        )


        gene_ratio = (
            k / N
        )


        background_ratio = (
            n / M
        )


        fold_enrichment = (
            gene_ratio
            / background_ratio
            if background_ratio > 0
            else np.nan
        )


        results.append(
            {
                "KEGG pathway": pathway,
                "Gene_count": k,
                "Query_size": N,
                "Pathway_size": n,
                "Background_size": M,
                "P_value": pvalue,
                "Fold_enrichment": fold_enrichment,
                "Genes": "; ".join(
                    sorted(
                        overlap
                    )
                ),
            }
        )


    if not results:

        return pd.DataFrame(
            columns=[
                "KEGG pathway",
                "Gene_count",
                "Query_size",
                "Pathway_size",
                "Background_size",
                "P_value",
                "Adjusted_P_value",
                "Fold_enrichment",
                "Genes",
            ]
        )


    df = pd.DataFrame(
        results
    )


    df[
        "Adjusted_P_value"
    ] = benjamini_hochberg(
        df[
            "P_value"
        ].values
    )


    df[
        "minus_log10_pvalue"
    ] = -np.log10(
        df[
            "P_value"
        ].clip(
            lower=1e-300
        )
    )


    df = df.sort_values(
        by=[
            "Adjusted_P_value",
            "P_value",
            "Gene_count"
        ],
        ascending=[
            True,
            True,
            False
        ]
    ).reset_index(
        drop=True
    )


    return df


# ============================================================
# Bubble plot
# ============================================================

def plot_enrichment(
    enrichment_df,
    title,
    output_file,
    top_n=10
):

    if enrichment_df.empty:

        print(
            f"WARNING: No pathways to plot for {title}"
        )

        return


    # --------------------------------------------------------
    # Use top pathways by p-value
    # --------------------------------------------------------

    plot_df = (
        enrichment_df
        .sort_values(
            by=[
                "P_value",
                "Gene_count"
            ],
            ascending=[
                True,
                False
            ]
        )
        .head(
            top_n
        )
        .copy()
    )


    # Reverse so most significant is at top
    plot_df = plot_df.iloc[
        ::-1
    ].reset_index(
        drop=True
    )


    # --------------------------------------------------------
    # Figure sizing
    # --------------------------------------------------------

    max_label_length = max(
        len(str(x))
        for x in plot_df[
            "KEGG pathway"
        ]
    )


    width = 8.5

    height = max(
        5.0,
        0.55
        * len(plot_df)
    )


    fig, ax = plt.subplots(
        figsize=(
            width,
            height
        )
    )


    y = np.arange(
        len(
            plot_df
        )
    )


    counts = plot_df[
        "Gene_count"
    ].values


    colors = plot_df[
        "minus_log10_pvalue"
    ].values


    # --------------------------------------------------------
    # Scale dot sizes
    # --------------------------------------------------------

    if counts.max() == counts.min():

        sizes = np.full(
            len(counts),
            180
        )

    else:

        sizes = (
            90
            + 420
            * (
                (
                    counts
                    - counts.min()
                )
                /
                (
                    counts.max()
                    - counts.min()
                )
            )
        )


    scatter = ax.scatter(
        counts,
        y,
        s=sizes,
        c=colors,
        cmap="coolwarm",
        edgecolors="none"
    )


    # --------------------------------------------------------
    # Y labels
    # --------------------------------------------------------

    ax.set_yticks(
        y
    )


    ax.set_yticklabels(
        plot_df[
            "KEGG pathway"
        ],
        fontsize=10
    )


    ax.set_xlabel(
        "Number of genes",
        fontsize=12
    )


    ax.set_title(
        title,
        fontsize=13
    )


    # --------------------------------------------------------
    # Grid
    # --------------------------------------------------------

    ax.grid(
        True,
        linewidth=0.6,
        alpha=0.3
    )


    ax.set_axisbelow(
        True
    )


    # --------------------------------------------------------
    # Color scale
    # --------------------------------------------------------

    cbar = fig.colorbar(
        scatter,
        ax=ax,
        fraction=0.045,
        pad=0.05
    )


    cbar.set_label(
        "−log10(p-value)",
        fontsize=10
    )


    # --------------------------------------------------------
    # Dot-size legend
    # --------------------------------------------------------

    unique_counts = sorted(
        set(
            counts
        )
    )


    if len(
        unique_counts
    ) > 3:

        legend_counts = [
            unique_counts[0],
            unique_counts[
                len(unique_counts) // 2
            ],
            unique_counts[-1]
        ]

    else:

        legend_counts = (
            unique_counts
        )


    handles = []


    for count in legend_counts:

        if counts.max() == counts.min():

            marker_size = math.sqrt(
                180
            )

        else:

            scaled_size = (
                90
                + 420
                * (
                    (
                        count
                        - counts.min()
                    )
                    /
                    (
                        counts.max()
                        - counts.min()
                    )
                )
            )

            marker_size = math.sqrt(
                scaled_size
            )


        handle = plt.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markersize=marker_size,
            markerfacecolor="black",
            markeredgecolor="black",
            label=str(count)
        )


        handles.append(
            handle
        )


    ax.legend(
        handles=handles,
        title="Genes",
        loc="center left",
        bbox_to_anchor=(
            1.18,
            0.25
        ),
        frameon=False
    )


    # --------------------------------------------------------
    # Leave room for pathway labels + legends
    # --------------------------------------------------------

    fig.subplots_adjust(
        left=0.43,
        right=0.80,
        top=0.90,
        bottom=0.12
    )


    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )


    png_file = output_file.with_suffix(
        ".png"
    )


    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )


    plt.close(
        fig
    )


    print(
        f"Saved plot: "
        f"{output_file}"
    )


# ============================================================
# Read cluster assignments
# ============================================================

def load_cluster_genes(
    cluster_file
):

    df = pd.read_csv(
        cluster_file
    )


    required = {
        "Geneid",
        "cluster"
    }


    missing = (
        required
        - set(
            df.columns
        )
    )


    if missing:

        raise ValueError(
            "Cluster file is missing: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )


    df[
        "cluster"
    ] = pd.to_numeric(
        df[
            "cluster"
        ],
        errors="raise"
    ).astype(
        int
    )


    df[
        "Geneid"
    ] = df[
        "Geneid"
    ].apply(
        normalize_gene_id
    )


    cluster1 = set(
        df.loc[
            df[
                "cluster"
            ] == 1,
            "Geneid"
        ]
    )


    cluster2 = set(
        df.loc[
            df[
                "cluster"
            ] == 2,
            "Geneid"
        ]
    )


    return (
        cluster1,
        cluster2
    )


# ============================================================
# Read significant genes from one timepoint
# ============================================================

def load_timepoint_genes(
    significant_file
):

    df = pd.read_csv(
        significant_file
    )


    possible_gene_cols = [
        "Geneid",
        "gene_id",
        "GeneID",
        "Unnamed: 0"
    ]


    gene_col = None


    for candidate in possible_gene_cols:

        if candidate in df.columns:

            gene_col = candidate

            break


    if gene_col is None:

        raise ValueError(
            "Could not find gene ID column in:\n"
            f"{significant_file}"
        )


    genes = set(
        df[
            gene_col
        ]
        .dropna()
        .apply(
            normalize_gene_id
        )
    )


    return genes


# ============================================================
# Run one analysis
# ============================================================

def analyze_gene_set(
    name,
    genes,
    background_genes,
    pathway_to_genes,
    output_dir,
    top_n
):

    print()
    print("========================================")
    print(name)
    print("========================================")


    genes_in_background = (
        genes
        & background_genes
    )


    print(
        f"Input genes: "
        f"{len(genes)}"
    )


    print(
        f"Genes in KEGG background: "
        f"{len(genes_in_background)}"
    )


    enrichment = run_kegg_enrichment(
        genes,
        background_genes,
        pathway_to_genes
    )


    safe_name = (
        name.lower()
        .replace(
            " ",
            "_"
        )
    )


    csv_file = (
        output_dir
        / f"{safe_name}_kegg_enrichment.csv"
    )


    enrichment.to_csv(
        csv_file,
        index=False
    )


    print(
        f"Saved table: "
        f"{csv_file}"
    )


    plot_enrichment(
        enrichment,
        title=f"{name} KEGG pathway enrichment",
        output_file=(
            output_dir
            / f"{safe_name}_kegg_enrichment.pdf"
        ),
        top_n=top_n
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "KEGG enrichment analysis for HELIOS NAD-RNA "
            "DTW clusters and a selected time point."
        )
    )


    parser.add_argument(
        "workflow_dir",
        help="Root tp_workflow directory."
    )


    parser.add_argument(
        "--cluster-dir",
        default=None,
        help=(
            "Step 16 clustering directory. Default: "
            "count_tables/intergenic/dtw_clustering_min3tp"
        )
    )


    parser.add_argument(
        "--timepoint",
        type=int,
        default=7,
        help=(
            "Time point for independent KEGG analysis. "
            "Default: 7."
        )
    )


    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help=(
            "Maximum pathways shown in each bubble plot. "
            "Default: 10."
        )
    )


    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help=(
            "Adjusted p-value threshold printed in summary. "
            "All pathways are still written to CSV. "
            "Default: 0.05."
        )
    )


    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. Default: "
            "<workflow>/kegg_enrichment"
        )
    )


    args = parser.parse_args()


    # ========================================================
    # Paths
    # ========================================================

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


    if args.cluster_dir:

        cluster_dir = Path(
            args.cluster_dir
        ).expanduser().resolve()

    else:

        cluster_dir = (
            intergenic_dir
            / "dtw_clustering_min3tp"
        )


    cluster_file = (
        cluster_dir
        / "nad_genes_clusters.csv"
    )


    tp_name = (
        f"tp{args.timepoint}"
    )


    tp_file = (
        intergenic_dir
        / tp_name
        / "significant_intergenic_Astart.csv"
    )


    gtf_file = (
        workflow_dir
        / "gtf"
        / "ncbi_dataset"
        / "data"
        / "GCF_000005845.2"
        / "genomic.gtf"
    )


    cache_dir = (
        workflow_dir
        / "kegg_cache"
    )


    if args.output_dir:

        output_dir = Path(
            args.output_dir
        ).expanduser().resolve()

    else:

        output_dir = (
            workflow_dir
            / "kegg_enrichment"
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # ========================================================
    # Validate
    # ========================================================

    for file_path in [
        cluster_file,
        tp_file,
        gtf_file
    ]:

        if not file_path.is_file():

            raise FileNotFoundError(
                f"Required file not found:\n"
                f"{file_path}"
            )


    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 17 - KEGG pathway enrichment")
    print("========================================")


    print(
        f"Cluster file:\n"
        f"  {cluster_file}"
    )


    print(
        f"Timepoint file:\n"
        f"  {tp_file}"
    )


    print(
        f"GTF:\n"
        f"  {gtf_file}"
    )


    print(
        f"Output:\n"
        f"  {output_dir}"
    )


    print("========================================")


    # ========================================================
    # KEGG
    # ========================================================

    pathway_names = (
        load_kegg_pathway_names(
            cache_dir
        )
    )


    (
        gene_to_pathways,
        pathway_to_genes
    ) = load_kegg_gene_pathways(
        cache_dir,
        pathway_names
    )


    # ========================================================
    # Background
    #
    # All GTF genes that also have a KEGG annotation.
    # ========================================================

    gtf_genes = load_gtf_genes(
        gtf_file
    )


    kegg_genes = set(
        gene_to_pathways.keys()
    )


    background_genes = (
        gtf_genes
        & kegg_genes
    )


    print(
        f"Genes in GTF: "
        f"{len(gtf_genes)}"
    )


    print(
        f"Genes with KEGG annotation: "
        f"{len(kegg_genes)}"
    )


    print(
        f"KEGG enrichment background: "
        f"{len(background_genes)} genes"
    )


    # ========================================================
    # Query sets
    # ========================================================

    (
        cluster1_genes,
        cluster2_genes
    ) = load_cluster_genes(
        cluster_file
    )


    tp_genes = load_timepoint_genes(
        tp_file
    )


    # ========================================================
    # Cluster 1
    # ========================================================

    analyze_gene_set(
        name="Cluster 1",
        genes=cluster1_genes,
        background_genes=background_genes,
        pathway_to_genes=pathway_to_genes,
        output_dir=output_dir,
        top_n=args.top_n
    )


    # ========================================================
    # Cluster 2
    # ========================================================

    analyze_gene_set(
        name="Cluster 2",
        genes=cluster2_genes,
        background_genes=background_genes,
        pathway_to_genes=pathway_to_genes,
        output_dir=output_dir,
        top_n=args.top_n
    )


    # ========================================================
    # tp7
    # ========================================================

    analyze_gene_set(
        name=tp_name,
        genes=tp_genes,
        background_genes=background_genes,
        pathway_to_genes=pathway_to_genes,
        output_dir=output_dir,
        top_n=args.top_n
    )


    print()
    print("========================================")
    print("Step 17 completed.")
    print("========================================")

    print(
        f"Cluster 1 genes: "
        f"{len(cluster1_genes)}"
    )

    print(
        f"Cluster 2 genes: "
        f"{len(cluster2_genes)}"
    )

    print(
        f"{tp_name} genes: "
        f"{len(tp_genes)}"
    )

    print()

    print(
        f"Output directory:\n"
        f"{output_dir}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
