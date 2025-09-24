#!/usr/bin/env python3
import argparse
import ast
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from tslearn.clustering import TimeSeriesKMeans
from tslearn.utils import to_time_series_dataset

# Golden ratio for figure aspect
PHI = (1 + np.sqrt(5)) / 2

# -- CONFIGURATION (defaults) --
DEFAULT_INFILE = Path(
    "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_Astart/"
    "common_nad_genes_across_timepoints_with_stats.csv"
)
DEFAULT_OUTDIR = Path(
    "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/"
    "table_Astart"
)
TP_ORDER = [f"tp{i}" for i in range(1, 17)]  # tp1..tp16

def parse_args():
    p = argparse.ArgumentParser(
        description="Cluster all genes (DTW, k=3) from TimePoints avg values; plot heatmap sorted by cluster."
    )
    p.add_argument("--infile", type=Path, default=DEFAULT_INFILE,
                   help="CSV with columns: Geneid, TimePoints, NumTimePoints")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                   help="Output directory for plots and matrices")
    p.add_argument("--n_clusters", type=int, default=3,
                   help="Number of DTW clusters (default=3)")
    return p.parse_args()

def parse_timepoints_cell(cell: str) -> dict:
    """
    Parse a cell like:
      "['tp1, 54455.75, 12651.32', 'tp2, 38291.20, 5665.01', ...]"
    Return dict { 'tp1': 54455.75, 'tp2': 38291.20, ... }
    If a value is missing or 'NA', store np.nan.
    """
    vals = {}
    if pd.isna(cell):
        return vals
    try:
        items = ast.literal_eval(str(cell))
    except Exception:
        return vals

    for item in items:
        parts = [p.strip() for p in str(item).split(",")]
        if not parts:
            continue
        tp = parts[0]
        # first numeric after tp is avg
        avg = np.nan
        for p in parts[1:]:
            try:
                avg = float(p)
                break
            except Exception:
                continue
        vals[tp] = avg
    return vals
  def build_matrix_from_common(file_path: Path) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    if not {"Geneid", "TimePoints"}.issubset(df.columns):
        raise ValueError("Input must contain 'Geneid' and 'TimePoints' columns")

    rows = []
    for _, r in df.iterrows():
        gene = str(r["Geneid"])
        tp_map = parse_timepoints_cell(r["TimePoints"])
        row = {"Geneid": gene}
        for tp in TP_ORDER:
            row[tp] = tp_map.get(tp, np.nan)
        rows.append(row)

    mat = pd.DataFrame(rows).set_index("Geneid")
    return mat  # raw avg values per tp (could contain NaNs)

def zscore_per_row_from_avg(mat: pd.DataFrame) -> pd.DataFrame:
    """
    Fill NaNs with 0, log10(avg+1), then z-score per gene (row).
    """
    mat_filled = mat.fillna(0.0)
    log_mat = np.log10(mat_filled + 1.0)
    mu = log_mat.mean(axis=1)
    sd = log_mat.std(axis=1).replace(0, 1)  # avoid divide-by-zero
    z = log_mat.sub(mu, axis=0).div(sd, axis=0)
    return z

def dtw_cluster(expr_z: pd.DataFrame, n_clusters: int, random_state: int = 0):
    """
    Cluster rows (genes) using TimeSeriesKMeans with DTW.
    """
    ts_data = to_time_series_dataset(expr_z.values)  # shape: (n_genes, n_timepoints, 1)
    model = TimeSeriesKMeans(
        n_clusters=n_clusters,
        metric="dtw",
        random_state=random_state,
        n_init=2,
    )
    labels = model.fit_predict(ts_data)  # 0..k-1
    return labels, model

def plot_heatmap(expr_z_sorted: pd.DataFrame, labels_sorted: np.ndarray, outdir: Path, title: str):
    outdir.mkdir(parents=True, exist_ok=True)
    # Map row colors by cluster (1..3)
    labels_1based = labels_sorted + 1
    color_map = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c"}
    row_colors = pd.Series(labels_1based, index=expr_z_sorted.index).map(color_map)

    n_genes = expr_z_sorted.shape[0]
    height = max(4, 0.18 * n_genes)  # scale with number of genes
    width = max(6, height * PHI)

    g = sns.clustermap(
        expr_z_sorted,
        cmap="viridis",
        col_cluster=False,
        row_cluster=False,
        yticklabels=False,
        row_colors=row_colors,
        cbar_kws={"label": "Z-score log10(avg + 1)", "shrink": 0.5},
        figsize=(width, height)
    )
    g.ax_heatmap.set_title(title)
    g.ax_heatmap.set_xlabel("Timepoint")
    g.ax_heatmap.set_ylabel("Genes (grouped by cluster)")
    out_pdf = outdir / "all_genes_heatmap_clusters.pdf"
    g.savefig(out_pdf, dpi=300)
    plt.close()
    print(f"Saved heatmap to: {out_pdf}")
  def main():
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Build gene x timepoint matrix from 'TimePoints' avg values
    mat_avg = build_matrix_from_common(args.infile)
    mat_avg.to_csv(outdir / "all_genes_avg_matrix.tsv", sep="\t")

    # 2) Z-score per gene after log10(avg+1)
    expr_z = zscore_per_row_from_avg(mat_avg)
    expr_z = expr_z[TP_ORDER]  # ensure ordered columns
    expr_z.to_csv(outdir / "all_genes_avg_matrix_zscore.tsv", sep="\t")

    # 3) DTW clustering
    labels, model = dtw_cluster(expr_z, n_clusters=args.n_clusters)

    # 4) Sort rows by cluster label and (optionally) by mean within cluster for neatness
    sort_index = np.lexsort((expr_z.mean(axis=1).values, labels))
    expr_z_sorted = expr_z.iloc[sort_index]
    labels_sorted = labels[sort_index]

    # 5) Save cluster assignments (1-based labels)
    cluster_df = pd.DataFrame({
        "Geneid": expr_z_sorted.index,
        "cluster": labels_sorted + 1
    })
    cluster_df.to_csv(outdir / "all_genes_clusters.csv", index=False)

    # 6) Cluster average trajectories
    avg_traj = expr_z.assign(cluster=labels).groupby("cluster").mean()
    avg_traj.index = avg_traj.index + 1  # 1-based for readability
    avg_traj.to_csv(outdir / "cluster_average_zscores.tsv", sep="\t")

    # Plot average trajectories vs time (minutes; 30 min spacing)
    time_minutes = [int(tp.lstrip("tp")) * 30 for tp in expr_z.columns]
    plt.figure(figsize=(6, 4))
    for c in avg_traj.index:
        plt.plot(time_minutes, avg_traj.loc[c], marker='o', label=f"Cluster {c}")
    plt.legend(fontsize='small', title='Cluster', title_fontsize='small')
    plt.xlabel("Time (min)")
    plt.ylabel("Average Z-score")
    plt.title(f"Average trajectories (k={args.n_clusters})")
    plt.xticks(time_minutes)
    plt.tight_layout()
    plt.savefig(outdir / "cluster_average_trajectories.pdf", dpi=300)
    plt.close()

    # 7) Heatmap (rows sorted by cluster)
    title = f"All genes (N={expr_z.shape[0]}) ?~@~T DTW clusters (k={args.n_clusters})"
    plot_heatmap(expr_z_sorted, labels_sorted, outdir, title)

if __name__ == "__main__":
    main()
