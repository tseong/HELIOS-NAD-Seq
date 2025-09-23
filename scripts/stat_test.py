#!/usr/bin/env python3
import os
import pickle as pkl
import pandas as pd

from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Set up inference engine
inference = DefaultInference(n_cpus=8)

# Base directory containing the tp1?~@~Stp16 subfolders
base_directory = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_Astart"

# Whether to save results to CSV
SAVE = True

# List of time-point folders
TP_LIST = [f"tp{i}" for i in range(1, 17)]

for tp in TP_LIST:
    tp_dir = os.path.join(base_directory, tp)
    if not os.path.isdir(tp_dir):
        print(f"[{tp}] Directory not found, skipping.")
        continue

    # Find the first *.pkl file in this subfolder
    dds_file = next((f for f in os.listdir(tp_dir) if f.endswith("_dds.pkl")), None)
    if not dds_file:
        print(f"[{tp}] No *_dds.pkl file found in {tp_dir}, skipping.")
        continue

    dds_path = os.path.join(tp_dir, dds_file)
    print(f"[{tp}] Loading DSS object from {dds_path}")
    with open(dds_path, "rb") as f:
        dds = pkl.load(f)

    # Run statistical testing
    print(f"[{tp}] Running DeseqStats")
    stat_res = DeseqStats(
        dds,
        alpha=0.05,
        cooks_filter=False,
        independent_filter=True
    )
    stat_res.run_wald_test()

    if stat_res.cooks_filter:
        stat_res._cooks_filtering()
    if stat_res.independent_filter:
        stat_res._independent_filtering()
    else:
        stat_res._p_value_adjustment()

    stat_res.summary()
    results_df = stat_res.results_df

    # Write out CSV
    output_name = dds_file.replace("_dds.pkl", "_results.csv")
    output_path = os.path.join(tp_dir, output_name)
    if SAVE:
        results_df.to_csv(output_path, index=True)
        print(f"[{tp}] Saved results to {output_path}")
