#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

from tslearn.clustering import TimeSeriesKMeans
from tslearn.utils import to_time_series_dataset
from tslearn.metrics import cdist_dtw

from sklearn.metrics import silhouette_score


# ============================================================
# HELIOS NAD-Seq pipeline
# Step 16: DTW clustering of common NAD genes
#
# Usage:
#
# python scripts/16.cluster_common_nad_genes.py \
#     <INTERGENIC_DIR>
#
# Example:
#
# python scripts/16.cluster_common_nad_genes.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow/count_tables/intergenic
#
# Optional:
#   --min-enriched-timepoints 3
#   --n-clusters 3
#   --k-min 2
#   --k-max 8
#   --start-tp 1
#   --end-tp 16
#   --random-state 0
#   --outdir PATH
# ============================================================


PHI = (1 + np.sqrt(5)) / 2


POSITIVE_BARCODES = [
    "bc01",
    "bc02",
    "bc03",
    "bc04",
]


CLUSTER_COLORS = [
    "#1f77b4",  # C1 blue
    "#ff7f0e",  # C2 orange
    "#2ca02c",  # C3 green
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
]


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Cluster HELIOS NAD genes using DTW and evaluate "
            "different cluster numbers using silhouette scores."
        )
    )

    parser.add_argument(
        "input_dir",
        help=(
            "Intergenic count directory containing "
            "common_nad_genes_across_timepoints.csv and "
            "tp1-tp16 folders."
        )
    )

    parser.add_argument(
        "--min-enriched-timepoints",
        type=int,
        default=3,
        help=(
            "Minimum number of enriched time points required. "
            "Default: 3."
        )
    )

    parser.add_argument(
        "--n-clusters",
        type=int,
        default=3,
        help="Final number of DTW clusters. Default: 3."
    )

    parser.add_argument(
        "--k-min",
        type=int,
        default=2,
        help="Minimum k for silhouette testing. Default: 2."
    )

    parser.add_argument(
        "--k-max",
        type=int,
        default=8,
        help="Maximum k for silhouette testing. Default: 8."
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
        "--random-state",
        type=int,
        default=0
    )

    parser.add_argument(
        "--outdir",
        default=None
    )

    return parser.parse_args()


# ============================================================
# Load common genes
# ============================================================

def load_common_genes(
    common_file,
    min_enriched_timepoints
):

    df = pd.read_csv(
        common_file
    )


    required = {
        "Geneid",
        "n_timepoints",
    }


    missing = (
        required
        - set(df.columns)
    )


    if missing:

        raise ValueError(
            "common_nad_genes_across_timepoints.csv is missing: "
            + ", ".join(sorted(missing))
        )


    df["n_timepoints"] = pd.to_numeric(
        df["n_timepoints"],
        errors="raise"
    )


    n_before = len(df)


    df = df[
        df["n_timepoints"]
        >= min_enriched_timepoints
    ].copy()


    print(
        f"Genes before filtering: {n_before}"
    )

    print(
        f"Genes enriched in >= "
        f"{min_enriched_timepoints} time points: "
        f"{len(df)}"
    )


    genes = (
        df["Geneid"]
        .dropna()
        .astype(str)
        .str.strip()
    )


    genes = sorted(
        set(
            gene
            for gene in genes
            if gene
        )
    )


    if not genes:

        raise ValueError(
            "No genes remain after filtering."
        )


    return genes


# ============================================================
# Construct gene x time matrix
# ============================================================

def build_average_matrix(
    input_dir,
    common_genes,
    tp_order
):

    common_gene_set = set(
        common_genes
    )


    matrix = pd.DataFrame(
        index=common_genes,
        columns=tp_order,
        dtype=float
    )


    for tp in tp_order:

        count_file = (
            input_dir
            / tp
            / "nad_genes_readCount.csv"
        )


        if not count_file.is_file():

            raise FileNotFoundError(
                f"Missing count file:\n"
                f"{count_file}"
            )


        print(
            f"Reading {tp}: {count_file}"
        )


        df = pd.read_csv(
            count_file
        )


        required = {
            "Geneid",
            *POSITIVE_BARCODES
        }


        missing = (
            required
            - set(df.columns)
        )


        if missing:

            raise ValueError(
                f"{count_file} is missing: "
                + ", ".join(sorted(missing))
            )


        df["Geneid"] = (
            df["Geneid"]
            .astype(str)
            .str.strip()
        )


        df = df[
            df["Geneid"].isin(
                common_gene_set
            )
        ].copy()


        df[
            POSITIVE_BARCODES
        ] = df[
            POSITIVE_BARCODES
        ].apply(
            pd.to_numeric,
            errors="raise"
        )


        # Mean of bc01-bc04
        df["mean_positive"] = (
            df[
                POSITIVE_BARCODES
            ]
            .mean(axis=1)
        )


        gene_means = (
            df.groupby("Geneid")[
                "mean_positive"
            ]
            .mean()
        )


        matrix.loc[
            gene_means.index,
            tp
        ] = gene_means


    return matrix


# ============================================================
# Transform and Z-score
# ============================================================

def transform_and_zscore(
    matrix
):

    matrix_filled = matrix.fillna(
        0.0
    )


    log_matrix = np.log10(
        matrix_filled + 1.0
    )


    row_mean = log_matrix.mean(
        axis=1
    )

    row_sd = log_matrix.std(
        axis=1
    )


    variable_mask = (
        row_sd > 0
    )


    removed = int(
        (~variable_mask).sum()
    )


    if removed > 0:

        print(
            f"Removing {removed} genes with "
            "zero variance across time."
        )


    log_matrix = log_matrix.loc[
        variable_mask
    ]


    row_mean = row_mean.loc[
        variable_mask
    ]

    row_sd = row_sd.loc[
        variable_mask
    ]


    z_matrix = (
        log_matrix
        .sub(
            row_mean,
            axis=0
        )
        .div(
            row_sd,
            axis=0
        )
    )


    return (
        log_matrix,
        z_matrix
    )


# ============================================================
# DTW clustering
# ============================================================

def run_dtw_clustering(
    z_matrix,
    n_clusters,
    random_state
):

    if len(z_matrix) <= n_clusters:

        raise ValueError(
            f"Only {len(z_matrix)} genes are available, "
            f"but n_clusters={n_clusters}."
        )


    ts_data = to_time_series_dataset(
        z_matrix.values
    )


    model = TimeSeriesKMeans(
        n_clusters=n_clusters,
        metric="dtw",
        random_state=random_state,
        n_init=5,
        max_iter=50
    )


    labels = model.fit_predict(
        ts_data
    )


    return (
        labels,
        model,
        ts_data
    )


# ============================================================
# Silhouette
# ============================================================

def calculate_silhouette(
    ts_data,
    labels,
    distance_matrix=None
):

    if len(
        np.unique(labels)
    ) < 2:

        return np.nan


    if distance_matrix is None:

        distance_matrix = cdist_dtw(
            ts_data
        )


    return silhouette_score(
        distance_matrix,
        labels,
        metric="precomputed"
    )


# ============================================================
# Test multiple k values
# ============================================================

def evaluate_k_values(
    z_matrix,
    k_min,
    k_max,
    random_state
):

    ts_data = to_time_series_dataset(
        z_matrix.values
    )


    print()
    print(
        "Calculating pairwise DTW distance matrix..."
    )


    distance_matrix = cdist_dtw(
        ts_data
    )


    results = []


    print()
    print("========================================")
    print("DTW silhouette scores by k")
    print("========================================")


    for k in range(
        k_min,
        k_max + 1
    ):

        if len(z_matrix) <= k:

            continue


        model = TimeSeriesKMeans(
            n_clusters=k,
            metric="dtw",
            random_state=random_state,
            n_init=5,
            max_iter=50
        )


        labels = model.fit_predict(
            ts_data
        )


        if len(
            np.unique(labels)
        ) < 2:

            score = np.nan

        else:

            score = silhouette_score(
                distance_matrix,
                labels,
                metric="precomputed"
            )


        results.append(
            {
                "k": k,
                "silhouette_score": score
            }
        )


        if np.isnan(score):

            print(
                f"k={k}: silhouette=NA"
            )

        else:

            print(
                f"k={k}: silhouette={score:.4f}"
            )


    print("========================================")


    return pd.DataFrame(
        results
    )


# ============================================================
# Silhouette plot
# ============================================================

def plot_silhouette_scores(
    scores_df,
    out_file
):

    valid = scores_df.dropna(
        subset=[
            "silhouette_score"
        ]
    )


    if valid.empty:

        return


    fig, ax = plt.subplots(
        figsize=(6, 4)
    )


    ax.plot(
        valid["k"],
        valid["silhouette_score"],
        marker="o"
    )


    ax.set_xlabel(
        "Number of clusters (k)"
    )

    ax.set_ylabel(
        "DTW silhouette score"
    )

    ax.set_title(
        "Cluster validation"
    )


    ax.set_xticks(
        valid["k"]
    )


    fig.tight_layout()


    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight"
    )


    plt.close(
        fig
    )


# ============================================================
# Heatmap
# ============================================================

def plot_heatmap(
    z_sorted,
    labels_sorted,
    out_pdf,
    min_enriched_timepoints,
    n_clusters
):

    n_genes = len(
        z_sorted
    )


    # --------------------------------------------------------
    # Landscape golden-ratio figure
    # --------------------------------------------------------

    width = 12.0
    height = width / PHI


    fig = plt.figure(
        figsize=(
            width,
            height
        )
    )


    # --------------------------------------------------------
    # Main layout:
    #
    # cluster labels | colored strip | heatmap
    #
    # Extra room on lower left is reserved for the
    # boxed Z-score legend.
    # --------------------------------------------------------

    grid = fig.add_gridspec(
        nrows=1,
        ncols=3,
        width_ratios=[
            0.75,
            0.28,
            10
        ],
        left=0.11,
        right=0.97,
        bottom=0.16,
        top=0.89,
        wspace=0.015
    )


    label_ax = fig.add_subplot(
        grid[0, 0]
    )

    cluster_ax = fig.add_subplot(
        grid[0, 1]
    )

    heat_ax = fig.add_subplot(
        grid[0, 2]
    )


    # ========================================================
    # Heatmap
    # ========================================================

    abs_max = float(
        np.nanmax(
            np.abs(
                z_sorted.values
            )
        )
    )


    # Use symmetric Z-score limits
    vmax = min(
        2.5,
        max(
            2.0,
            abs_max
        )
    )

    vmin = -vmax


    image = heat_ax.imshow(
        z_sorted.values,
        aspect="auto",
        interpolation="nearest",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax
    )


    # --------------------------------------------------------
    # Time axis
    #
    # tp2 = 1 h
    # tp4 = 2 h
    # ...
    # tp16 = 8 h
    # --------------------------------------------------------

    tick_positions = []
    tick_labels = []


    for i, tp in enumerate(
        z_sorted.columns
    ):

        tp_number = int(
            tp.replace(
                "tp",
                ""
            )
        )


        hours = (
            tp_number
            * 0.5
        )


        if hours.is_integer():

            tick_positions.append(
                i
            )

            tick_labels.append(
                str(
                    int(hours)
                )
            )


    heat_ax.set_xticks(
        tick_positions
    )

    heat_ax.set_xticklabels(
        tick_labels,
        fontsize=11
    )


    heat_ax.set_yticks(
        []
    )


    heat_ax.set_xlabel(
        "Time (h)",
        fontsize=13
    )


    heat_ax.set_title(
        "Heat map showing dynamics of NAD-capped RNA expression",
        fontsize=14,
        pad=12
    )


    # Remove unnecessary spines
    heat_ax.spines[
        "top"
    ].set_visible(
        False
    )

    heat_ax.spines[
        "right"
    ].set_visible(
        False
    )


    # ========================================================
    # Cluster color strip
    # ========================================================

    labels_1based = (
        labels_sorted
        + 1
    )


    cluster_cmap = ListedColormap(
        CLUSTER_COLORS[
            :n_clusters
        ]
    )


    cluster_values = (
        labels_1based
        - 1
    ).reshape(
        -1,
        1
    )


    cluster_ax.imshow(
        cluster_values,
        aspect="auto",
        interpolation="nearest",
        cmap=cluster_cmap,
        vmin=0,
        vmax=max(
            n_clusters - 1,
            1
        )
    )


    cluster_ax.set_xticks(
        []
    )

    cluster_ax.set_yticks(
        []
    )


    for spine in (
        cluster_ax.spines.values()
    ):

        spine.set_visible(
            False
        )


    # ========================================================
    # C1 / C2 / C3 labels
    # ========================================================

    label_ax.set_xlim(
        0,
        1
    )

    label_ax.set_ylim(
        n_genes - 0.5,
        -0.5
    )

    label_ax.axis(
        "off"
    )


    unique_clusters = np.unique(
        labels_1based
    )


    for cluster in unique_clusters:

        positions = np.where(
            labels_1based == cluster
        )[0]


        if len(positions) == 0:

            continue


        start = positions[0]
        end = positions[-1]

        center = (
            start + end
        ) / 2


        label_ax.text(
            0.88,
            center,
            f"C{cluster}",
            ha="right",
            va="center",
            fontsize=13
        )


        # Cluster boundary
        if end < n_genes - 1:

            heat_ax.axhline(
                end + 0.5,
                linewidth=0.8,
                color="white"
            )

            cluster_ax.axhline(
                end + 0.5,
                linewidth=0.8,
                color="white"
            )


    # ========================================================
    # BOXED Z-SCORE LEGEND
    #
    # This is deliberately created as its own axes rather
    # than as an inset of the cluster-label axis.
    #
    # Everything is therefore contained inside one box:
    #
    #   Z-score
    #   log10(avg + 1)
    #   gradient
    #   2
    #   1
    #   0
    #  -1
    #  -2
    #
    # ========================================================

    legend_ax = fig.add_axes(
        [
            0.018,   # left
            0.055,   # bottom
            0.075,   # width
            0.22     # height
        ]
    )


    # White background with visible gray box
    legend_ax.set_facecolor(
        "white"
    )


    for spine in (
        legend_ax.spines.values()
    ):

        spine.set_visible(
            True
        )

        spine.set_linewidth(
            0.8
        )

        spine.set_edgecolor(
            "0.55"
        )


    legend_ax.set_xticks(
        []
    )

    legend_ax.set_yticks(
        []
    )


    legend_ax.set_xlim(
        0,
        1
    )

    legend_ax.set_ylim(
        0,
        1
    )


    # --------------------------------------------------------
    # Legend title
    # --------------------------------------------------------

    legend_ax.text(
        0.5,
        0.93,
        "Z-score",
        ha="center",
        va="top",
        fontsize=8
    )


    legend_ax.text(
        0.5,
        0.82,
        "log10(avg + 1)",
        ha="center",
        va="top",
        fontsize=6.5
    )


    # --------------------------------------------------------
    # Gradient area inside box
    # --------------------------------------------------------

    gradient_ax = legend_ax.inset_axes(
        [
            0.28,
            0.12,
            0.28,
            0.58
        ]
    )


    gradient = np.linspace(
        vmin,
        vmax,
        256
    ).reshape(
        -1,
        1
    )


    gradient_ax.imshow(
        gradient,
        aspect="auto",
        origin="lower",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax
    )


    gradient_ax.set_xticks(
        []
    )


    # --------------------------------------------------------
    # Z-score ticks inside legend box
    # --------------------------------------------------------

    desired_ticks = [
        -2,
        -1,
        0,
        1,
        2
    ]


    ticks = [
        tick
        for tick in desired_ticks
        if (
            tick >= vmin
            and tick <= vmax
        )
    ]


    tick_positions_gradient = [
        (
            tick - vmin
        )
        / (
            vmax - vmin
        )
        * 255

        for tick in ticks
    ]


    gradient_ax.set_yticks(
        tick_positions_gradient
    )


    gradient_ax.set_yticklabels(
        [
            str(tick)
            for tick in ticks
        ],
        fontsize=7
    )


    gradient_ax.yaxis.tick_right()


    gradient_ax.tick_params(
        axis="y",
        length=2,
        pad=2
    )


    for spine in (
        gradient_ax.spines.values()
    ):

        spine.set_visible(
            False
        )


    # ========================================================
    # Save
    # ========================================================

    fig.savefig(
        out_pdf,
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )


    out_png = out_pdf.with_suffix(
        ".png"
    )


    fig.savefig(
        out_png,
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )


    plt.close(
        fig
    )


    print(
        f"Saved heatmap PDF: "
        f"{out_pdf}"
    )

    print(
        f"Saved heatmap PNG: "
        f"{out_png}"
    )


# ============================================================
# Average trajectories
# ============================================================

def plot_cluster_trajectories(
    average_trajectories,
    out_file
):

    time_hours = [
        int(
            tp.replace(
                "tp",
                ""
            )
        ) * 0.5

        for tp
        in average_trajectories.columns
    ]


    fig, ax = plt.subplots(
        figsize=(8, 5)
    )


    for cluster in (
        average_trajectories.index
    ):

        ax.plot(
            time_hours,
            average_trajectories.loc[
                cluster
            ].values,
            marker="o",
            label=f"Cluster {cluster}"
        )


    ax.set_xlabel(
        "Time (h)"
    )

    ax.set_ylabel(
        "Average Z-score"
    )

    ax.set_title(
        "Average NAD-gene trajectories"
    )


    ax.legend(
        title="Cluster"
    )


    fig.tight_layout()


    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight"
    )


    plt.close(
        fig
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()


    input_dir = Path(
        args.input_dir
    ).expanduser().resolve()


    if not input_dir.is_dir():

        raise FileNotFoundError(
            f"Input directory does not exist:\n"
            f"{input_dir}"
        )


    if args.min_enriched_timepoints < 1:

        raise ValueError(
            "--min-enriched-timepoints must be >= 1."
        )


    if args.n_clusters < 2:

        raise ValueError(
            "--n-clusters must be >= 2."
        )


    if args.n_clusters > len(
        CLUSTER_COLORS
    ):

        raise ValueError(
            "Too many clusters for available colors."
        )


    if args.k_min < 2:

        raise ValueError(
            "--k-min must be >= 2."
        )


    if args.k_max < args.k_min:

        raise ValueError(
            "--k-max must be >= --k-min."
        )


    if args.start_tp > args.end_tp:

        raise ValueError(
            "--start-tp must be <= --end-tp."
        )


    tp_order = [
        f"tp{i}"
        for i in range(
            args.start_tp,
            args.end_tp + 1
        )
    ]


    common_file = (
        input_dir
        / "common_nad_genes_across_timepoints.csv"
    )


    if not common_file.is_file():

        raise FileNotFoundError(
            "Common NAD-gene file not found:\n"
            f"{common_file}"
        )


    if args.outdir:

        outdir = Path(
            args.outdir
        ).expanduser().resolve()

    else:

        outdir = (
            input_dir
            / (
                "dtw_clustering_min"
                f"{args.min_enriched_timepoints}tp"
            )
        )


    outdir.mkdir(
        parents=True,
        exist_ok=True
    )


    print("========================================")
    print("HELIOS NAD-Seq")
    print("Step 16 - DTW clustering")
    print("========================================")

    print(
        f"Input:              {input_dir}"
    )

    print(
        f"Minimum enrichment: "
        f"{args.min_enriched_timepoints} time points"
    )

    print(
        f"Final clustering k: "
        f"{args.n_clusters}"
    )

    print(
        f"Silhouette k range: "
        f"{args.k_min}-{args.k_max}"
    )

    print(
        f"Output:             {outdir}"
    )

    print("========================================")


    # ========================================================
    # 1. Genes
    # ========================================================

    common_genes = load_common_genes(
        common_file,
        args.min_enriched_timepoints
    )


    # ========================================================
    # 2. Count matrix
    # ========================================================

    average_matrix = build_average_matrix(
        input_dir,
        common_genes,
        tp_order
    )


    average_matrix.to_csv(
        outdir
        / "nad_genes_average_counts.tsv",
        sep="\t"
    )


    # ========================================================
    # 3. Transformation
    # ========================================================

    (
        log_matrix,
        z_matrix
    ) = transform_and_zscore(
        average_matrix
    )


    log_matrix.to_csv(
        outdir
        / "nad_genes_log10_counts.tsv",
        sep="\t"
    )


    z_matrix.to_csv(
        outdir
        / "nad_genes_zscores.tsv",
        sep="\t"
    )


    print(
        f"Genes available for clustering: "
        f"{len(z_matrix)}"
    )


    # ========================================================
    # 4. k evaluation
    # ========================================================

    k_scores_df = evaluate_k_values(
        z_matrix,
        args.k_min,
        args.k_max,
        args.random_state
    )


    k_scores_file = (
        outdir
        / "silhouette_scores_by_k.csv"
    )


    k_scores_df.to_csv(
        k_scores_file,
        index=False
    )


    plot_silhouette_scores(
        k_scores_df,
        outdir
        / "silhouette_scores_by_k.pdf"
    )


    valid_scores = (
        k_scores_df
        .dropna(
            subset=[
                "silhouette_score"
            ]
        )
    )


    if not valid_scores.empty:

        best_row = valid_scores.loc[
            valid_scores[
                "silhouette_score"
            ].idxmax()
        ]


        best_k = int(
            best_row["k"]
        )

        best_score = float(
            best_row[
                "silhouette_score"
            ]
        )


        print()
        print(
            f"Best silhouette score: "
            f"k={best_k}, "
            f"score={best_score:.4f}"
        )

    else:

        best_k = None
        best_score = np.nan


    # ========================================================
    # 5. Final clustering
    # ========================================================

    labels, model, ts_data = (
        run_dtw_clustering(
            z_matrix,
            args.n_clusters,
            args.random_state
        )
    )


    final_distance_matrix = cdist_dtw(
        ts_data
    )


    final_silhouette = calculate_silhouette(
        ts_data,
        labels,
        distance_matrix=final_distance_matrix
    )


    print()
    print(
        f"Final clustering silhouette "
        f"(k={args.n_clusters}): "
        f"{final_silhouette:.4f}"
    )


    # ========================================================
    # 6. Metrics
    # ========================================================

    metrics_df = pd.DataFrame(
        {
            "metric": [
                "final_dtw_silhouette",
                "best_tested_dtw_silhouette"
            ],

            "value": [
                final_silhouette,
                best_score
            ],

            "k": [
                args.n_clusters,
                best_k
            ],

            "n_genes": [
                len(z_matrix),
                len(z_matrix)
            ]
        }
    )


    metrics_df.to_csv(
        outdir
        / "clustering_metrics.csv",
        index=False
    )


    # ========================================================
    # 7. Sort strictly by cluster
    # ========================================================

    sort_index = np.argsort(
        labels,
        kind="stable"
    )


    z_sorted = z_matrix.iloc[
        sort_index
    ]


    labels_sorted = labels[
        sort_index
    ]


    # ========================================================
    # 8. Cluster assignments
    # ========================================================

    cluster_df = pd.DataFrame(
        {
            "Geneid": (
                z_sorted.index
            ),

            "cluster": (
                labels_sorted + 1
            )
        }
    )


    cluster_df.to_csv(
        outdir
        / "nad_genes_clusters.csv",
        index=False
    )


    # ========================================================
    # 9. Average trajectories
    # ========================================================

    trajectory_df = z_matrix.copy()


    trajectory_df[
        "cluster"
    ] = (
        labels + 1
    )


    average_trajectories = (
        trajectory_df
        .groupby(
            "cluster"
        )
        .mean()
    )


    average_trajectories.to_csv(
        outdir
        / "cluster_average_zscores.tsv",
        sep="\t"
    )


    # ========================================================
    # 10. Cluster sizes
    # ========================================================

    cluster_sizes = (
        cluster_df[
            "cluster"
        ]
        .value_counts()
        .sort_index()
    )


    print()
    print("Final cluster sizes:")


    for cluster, n_genes in (
        cluster_sizes.items()
    ):

        print(
            f"  Cluster {cluster}: "
            f"{n_genes} genes"
        )


    # ========================================================
    # 11. Trajectory plot
    # ========================================================

    plot_cluster_trajectories(
        average_trajectories,
        outdir
        / "cluster_average_trajectories.pdf"
    )


    # ========================================================
    # 12. Heatmap
    # ========================================================

    plot_heatmap(
        z_sorted,
        labels_sorted,
        outdir
        / "nad_genes_heatmap_clusters.pdf",
        args.min_enriched_timepoints,
        args.n_clusters
    )


    print()
    print("========================================")
    print("Step 16 completed.")
    print("----------------------------------------")

    print(
        f"Genes clustered: "
        f"{len(z_matrix)}"
    )

    print(
        f"Final k: "
        f"{args.n_clusters}"
    )

    print(
        f"Final silhouette: "
        f"{final_silhouette:.4f}"
    )


    if best_k is not None:

        print(
            f"Best tested k: "
            f"{best_k}"
        )

        print(
            f"Best tested silhouette: "
            f"{best_score:.4f}"
        )


    print()
    print(
        f"Output directory:\n"
        f"{outdir}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
