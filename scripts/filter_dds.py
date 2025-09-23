#!/usr/bin/env python3
import os
import pickle as pkl
import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

# Base paths
BASE_DIR    = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2"
TABLE_DIR   = os.path.join(BASE_DIR, "table_Astart")

# Loop through tp1?~@~Stp16
for tp in range(1, 17):
    subdir   = os.path.join(TABLE_DIR, f"tp{tp}")
    inp_csv  = os.path.join(subdir, "merged_by_barcode_Astart_readCount.csv")
    out_pkl  = os.path.join(subdir, "merged_by_barcode_Astart_readCount_dds.pkl")

    # Skip if input is missing
    if not os.path.exists(inp_csv):
        print(f"[tp{tp}] skipping?~@~Tno file at {inp_csv}")
        continue

    print(f"[tp{tp}] processing {inp_csv}")

    # --- Read & prepare counts matrix ---
    df = pd.read_csv(inp_csv, index_col=0)
    if df.shape[1] != 8:
        raise ValueError(f"[tp{tp}] expected 8 samples, found {df.shape[1]} columns")

    # Build metadata
    metadata = pd.DataFrame({
        "sample_id": df.columns,
        "conditions": ["Treated"] * 4 + ["Control"] * 4
    }, index=df.columns)

    # Numeric counts and filter low-expression
    counts = df.apply(pd.to_numeric, errors="coerce")
    keep   = counts.sum(axis=1) >= 10
    counts = counts.loc[keep].T  # transpose: samples?~Wgenes

    # --- DESeq2 setup & fit ---
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors="conditions",
        refit_cooks=False,
        inference=inference
    )
    dds.fit_size_factors(fit_type="ratio")
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    if dds.refit_cooks:
        dds.varm["replaced"] = np.zeros_like(dds.var_names, dtype=bool)
        dds.refit()

    # --- Save the DDS object ---
    with open(out_pkl, "wb") as f:
        pkl.dump(dds, f)
    print(f"[tp{tp}] saved DDS to {out_pkl}")
