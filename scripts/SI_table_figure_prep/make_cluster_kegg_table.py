#!/usr/bin/env python3

import argparse
import re
import urllib.request
from pathlib import Path

import pandas as pd


# ============================================================
# HELIOS NAD-Seq pipeline
# SI table preparation
#
# Generate a 4-column table:
#
#   Gene name
#   Gene biotype
#   KEGG pathway
#   Cluster
#
# Reports ALL KEGG pathways associated with each gene.
#
# Usage:
#
# python scripts/SI_table_figure_prep/make_cluster_kegg_table.py \
#     <WORKFLOW_DIR>
#
# Example:
#
# python scripts/SI_table_figure_prep/make_cluster_kegg_table.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow
#
# Optional:
#   --min-timepoints 3
#   --cluster-dir PATH
#   --output PATH
# ============================================================


# ------------------------------------------------------------
# GTF attribute parser
# ------------------------------------------------------------

ATTR_RE = re.compile(
    r'(\S+)\s+"([^"]+)"'
)


def parse_gtf_attributes(text):

    return {
        key: value
        for key, value
        in ATTR_RE.findall(text)
    }


# ------------------------------------------------------------
# Normalize gene identifier
# ------------------------------------------------------------

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
# Clean KEGG pathway name
# ============================================================

def clean_pathway_name(pathway_name):
    """
    Remove organism suffix from KEGG pathway names.

    Example:

        Pyruvate metabolism - Escherichia coli K-12 MG1655

    becomes:

        Pyruvate metabolism
    """

    pathway_name = str(
        pathway_name
    ).strip()


    pathway_name = re.sub(
        r"\s*-\s*Escherichia coli K-12 MG1655\s*$",
        "",
        pathway_name
    )


    return pathway_name.strip()


# ============================================================
# Load GTF annotations
# ============================================================

def load_gtf_annotations(gtf_file):

    annotation = {}


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


            if not gene_id and not locus_tag:
                continue


            normalized_ids = set()


            if gene_id:

                normalized_ids.add(
                    normalize_gene_id(
                        gene_id
                    )
                )


            if locus_tag:

                normalized_ids.add(
                    normalize_gene_id(
                        locus_tag
                    )
                )


            gene_name = (
                attrs.get("gene")
                or attrs.get("gene_name")
                or locus_tag
                or gene_id
                or "unknown"
            )


            gene_biotype = (
                attrs.get("gene_biotype")
                or attrs.get("gene_type")
                or attrs.get("transcript_biotype")
                or attrs.get("transcript_type")
                or "unknown"
            )


            info = {
                "gene_name": gene_name,
                "gene_biotype": gene_biotype,
                "locus_tag": (
                    normalize_gene_id(
                        locus_tag
                    )
                    if locus_tag
                    else None
                ),
            }


            for identifier in normalized_ids:

                annotation[
                    identifier
                ] = info


    return annotation


# ============================================================
# Download/cache KEGG text
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
        "Downloading KEGG annotation:"
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
            "Could not retrieve KEGG annotation and no "
            f"cache was available:\n{cache_file}\n\n"
            f"Original error: {exc}"
        )


    cache_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    cache_file.write_text(
        text
    )


    print(
        f"Cached KEGG data: "
        f"{cache_file}"
    )


    return text


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
            .strip()
            .replace(
                "path:",
                ""
            )
        )


        pathway_name = clean_pathway_name(
            parts[1]
        )


        pathway_names[
            pathway_id
        ] = pathway_name


    return pathway_names


# ============================================================
# Load KEGG gene -> pathway mapping
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


    for line in text.splitlines():

        if not line.strip():
            continue


        parts = line.split(
            "\t"
        )


        if len(parts) != 2:
            continue


        gene_id = normalize_gene_id(
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


        pathway_name = clean_pathway_name(
            pathway_name
        )


        gene_to_pathways.setdefault(
            gene_id,
            set()
        ).add(
            pathway_name
        )


    return gene_to_pathways


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prepare SI table containing gene name, "
            "gene biotype, all KEGG pathways and DTW cluster."
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
            "Minimum number of enriched time points. "
            "Default: 3."
        )
    )


    parser.add_argument(
        "--cluster-dir",
        default=None,
        help=(
            "Step 16 clustering directory."
        )
    )


    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output CSV. Default: "
            "<workflow>/SI_table_prep/"
            "nad_genes_kegg_clusters.csv"
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Paths
    # --------------------------------------------------------

    workflow_dir = Path(
        args.workflow_dir
    ).expanduser().resolve()


    if not workflow_dir.is_dir():

        raise FileNotFoundError(
            f"Workflow directory not found:\n"
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


    if args.cluster_dir:

        cluster_dir = Path(
            args.cluster_dir
        ).expanduser().resolve()

    else:

        cluster_dir = (
            intergenic_dir
            / (
                "dtw_clustering_min"
                f"{args.min_timepoints}tp"
            )
        )


    cluster_file = (
        cluster_dir
        / "nad_genes_clusters.csv"
    )


    cache_dir = (
        workflow_dir
        / "SI_table_prep"
        / "kegg_cache"
    )


    if args.output:

        output_file = Path(
            args.output
        ).expanduser().resolve()

    else:

        output_file = (
            workflow_dir
            / "SI_table_prep"
            / "nad_genes_kegg_clusters.csv"
        )


    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Validate
    # --------------------------------------------------------

    for path in [
        common_file,
        gtf_file,
        cluster_file,
    ]:

        if not path.is_file():

            raise FileNotFoundError(
                f"Required file not found:\n"
                f"{path}"
            )


    print("========================================")
    print("HELIOS NAD-Seq")
    print("SI table: KEGG + DTW clusters")
    print("========================================")

    print(
        f"Minimum time points: "
        f"{args.min_timepoints}"
    )

    print(
        f"Common genes:\n"
        f"  {common_file}"
    )

    print(
        f"Cluster assignments:\n"
        f"  {cluster_file}"
    )

    print(
        f"GTF:\n"
        f"  {gtf_file}"
    )

    print("========================================")


    # ========================================================
    # Load common genes
    # ========================================================

    common_df = pd.read_csv(
        common_file
    )


    required_common = {
        "Geneid",
        "n_timepoints"
    }


    missing = (
        required_common
        - set(
            common_df.columns
        )
    )


    if missing:

        raise ValueError(
            "Common-gene file is missing: "
            + ", ".join(
                sorted(missing)
            )
        )


    common_df[
        "n_timepoints"
    ] = pd.to_numeric(
        common_df[
            "n_timepoints"
        ],
        errors="raise"
    )


    common_df = common_df[
        common_df[
            "n_timepoints"
        ] >= args.min_timepoints
    ].copy()


    common_df[
        "normalized_gene_id"
    ] = common_df[
        "Geneid"
    ].apply(
        normalize_gene_id
    )


    eligible_genes = set(
        common_df[
            "normalized_gene_id"
        ]
    )


    print(
        f"Genes with n_timepoints >= "
        f"{args.min_timepoints}: "
        f"{len(eligible_genes)}"
    )


    # ========================================================
    # Load cluster assignments
    # ========================================================

    cluster_df = pd.read_csv(
        cluster_file
    )


    required_cluster = {
        "Geneid",
        "cluster"
    }


    missing = (
        required_cluster
        - set(
            cluster_df.columns
        )
    )


    if missing:

        raise ValueError(
            "Cluster file is missing: "
            + ", ".join(
                sorted(missing)
            )
        )


    cluster_df[
        "cluster"
    ] = pd.to_numeric(
        cluster_df[
            "cluster"
        ],
        errors="raise"
    ).astype(
        int
    )


    cluster_df[
        "normalized_gene_id"
    ] = cluster_df[
        "Geneid"
    ].apply(
        normalize_gene_id
    )


    cluster_df = cluster_df[
        cluster_df[
            "normalized_gene_id"
        ].isin(
            eligible_genes
        )
    ].copy()


    # --------------------------------------------------------
    # Only clusters 1 and 2
    # --------------------------------------------------------

    observed_clusters = set(
        cluster_df[
            "cluster"
        ].unique()
    )


    invalid_clusters = (
        observed_clusters
        - {1, 2}
    )


    if invalid_clusters:

        raise ValueError(
            "Cluster file contains clusters other than "
            "1 or 2. Run Step 16 using --n-clusters 2."
        )


    print(
        f"Genes in cluster file: "
        f"{len(cluster_df)}"
    )


    # ========================================================
    # GTF
    # ========================================================

    annotations = load_gtf_annotations(
        gtf_file
    )


    print(
        f"GTF annotations loaded: "
        f"{len(annotations)}"
    )


    # ========================================================
    # KEGG
    # ========================================================

    pathway_names = (
        load_kegg_pathway_names(
            cache_dir
        )
    )


    gene_to_pathways = (
        load_kegg_gene_pathways(
            cache_dir,
            pathway_names
        )
    )


    print(
        f"Genes with KEGG pathways: "
        f"{len(gene_to_pathways)}"
    )


    # ========================================================
    # Build output
    # ========================================================

    output_rows = []

    no_gtf = 0
    no_kegg = 0


    for _, row in cluster_df.iterrows():

        gene_id = row[
            "normalized_gene_id"
        ]


        cluster = int(
            row[
                "cluster"
            ]
        )


        meta = annotations.get(
            gene_id
        )


        if meta is None:

            gene_name = gene_id
            gene_biotype = "unknown"

            no_gtf += 1

        else:

            gene_name = meta[
                "gene_name"
            ]

            gene_biotype = meta[
                "gene_biotype"
            ]


        # ----------------------------------------------------
        # Determine KEGG ID
        # ----------------------------------------------------

        kegg_id = gene_id


        if (
            meta is not None
            and meta.get(
                "locus_tag"
            )
        ):

            kegg_id = meta[
                "locus_tag"
            ]


        pathways = gene_to_pathways.get(
            kegg_id,
            set()
        )


        # ----------------------------------------------------
        # Report ALL pathways
        # ----------------------------------------------------

        if pathways:

            cleaned_pathways = sorted(
                {
                    clean_pathway_name(
                        pathway
                    )
                    for pathway in pathways
                    if clean_pathway_name(
                        pathway
                    )
                }
            )


            pathway_text = "; ".join(
                cleaned_pathways
            )

        else:

            pathway_text = "NA"

            no_kegg += 1


        output_rows.append(
            {
                "Gene name": gene_name,
                "Gene biotype": gene_biotype,
                "KEGG pathway": pathway_text,
                "Cluster": cluster,
            }
        )


    # ========================================================
    # Final table
    # ========================================================

    output_df = pd.DataFrame(
        output_rows,
        columns=[
            "Gene name",
            "Gene biotype",
            "KEGG pathway",
            "Cluster",
        ]
    )


    output_df = (
        output_df
        .sort_values(
            by=[
                "Cluster",
                "Gene name"
            ],
            ascending=[
                True,
                True
            ]
        )
        .reset_index(
            drop=True
        )
    )


    # ========================================================
    # Write output
    # ========================================================

    output_df.to_csv(
        output_file,
        index=False
    )


    # ========================================================
    # Summary
    # ========================================================

    print()
    print("========================================")
    print("SI table completed")
    print("----------------------------------------")

    print(
        f"Genes written: "
        f"{len(output_df)}"
    )


    for cluster in [
        1,
        2
    ]:

        n = int(
            (
                output_df[
                    "Cluster"
                ] == cluster
            ).sum()
        )

        print(
            f"Cluster {cluster}: "
            f"{n} genes"
        )


    print(
        f"Genes without GTF annotation: "
        f"{no_gtf}"
    )

    print(
        f"Genes without KEGG pathway: "
        f"{no_kegg}"
    )


    print()
    print(
        f"Output:\n"
        f"{output_file}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
