#!/usr/bin/env python3

import os
import pandas as pd

# === Configuration ===
# Path to the main TSV/space?~@~Pdelimited file you want to normalize
MAIN_FILE = '/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_filtered/nad_genes_normalized/nad_genes_across_tp_with_readCount_helios.tsv'

# Normalization files (as provided)
NORM_FILES = {
    'read_depth': '/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/trimmed_trimmomatic/'
                  'read_depth_by_timepoint_all_reads.csv',
    'assigned_reads': '/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/'
                      'assigned_reads.csv'
}

# Sampling summary file (10 rows: Sampling_1 ?~@? Sampling_10; columns tp1?~@?tp16, SUM)
SAMPLING_FILE = '/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/' \
                'common_genes_tp_least10_assigned_read_normalized/' \
                'tp_log2FC_common_genes_assigned_reads_normalized_summary.csv'

# Output folder for all normalized results
OUTPUT_DIR = '/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_filtered'
# =====================

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Columns we expect in the main file
EXPECTED_MAIN_COLS = {'gene_name', 'timepoint', 'gene_biotype',
                      '3PAB_rep1', '3PAB_rep2', '3PAB_rep3', '3PAB_rep4'}

# Function to load the main file, trying tab first then whitespace
def load_main_dataframe(path):
    # Attempt to read as tab-separated
    df = pd.read_csv(path, sep='\t')
    if EXPECTED_MAIN_COLS.issubset(df.columns):
        return df

    # If columns are missing, try reading as whitespace?~@~Pseparated
    df = pd.read_csv(path, sep=r'\s+', engine='python')
    if EXPECTED_MAIN_COLS.issubset(df.columns):
        return df

    # If still missing, report exactly which columns are missing
    missing = EXPECTED_MAIN_COLS - set(df.columns)
    raise ValueError(f"Missing required columns in {path}: {missing}")

# Load the main file into a DataFrame
df = load_main_dataframe(MAIN_FILE)

# We'll keep each initial normalization in memory so we can apply
# the 10 sampling-based normalizations later.
normalized_dfs = {}

# List of the 3PAB columns to normalize
REPS = ['3PAB_rep1', '3PAB_rep2', '3PAB_rep3', '3PAB_rep4']
# Step 1: Apply the two ?~@~\read_depth?~@~] and ?~@~\assigned_reads?~@~] normalizations
for norm_label, norm_path in NORM_FILES.items():
    # Load normalization CSV (assumed to have columns tp1 ?~@? tp16)
    norm_df = pd.read_csv(norm_path)
    if norm_df.shape[0] < 1:
        raise ValueError(f"Normalization file {norm_path} is empty or has no rows.")
    norm_row = norm_df.iloc[0]  # take the first row

    if 'tp1' not in norm_row.index:
        raise KeyError(f"'tp1' column not found in normalization file {norm_path}")

    # Compute normalization factors: for each tpX, factor = tpX / tp1
    factors = {}
    for tp_col in norm_row.index:
        if not tp_col.startswith('tp'):
            continue
        try:
            numerator = float(norm_row[tp_col])
            denominator = float(norm_row['tp1'])
            if denominator == 0:
                raise ZeroDivisionError(f"tp1 is zero in {norm_path}, cannot divide.")
            factors[tp_col] = numerator / denominator
        except ValueError:
            raise ValueError(f"Non numeric value in column '{tp_col}' of {norm_path}")

    # Helper: map a timepoint value to its normalization factor
    def get_factor(tp):
        tp_str = str(tp).strip()
        if not tp_str.startswith('tp'):
            tp_str = 'tp' + tp_str
        return factors.get(tp_str, 1.0)

    # Copy and normalize
    df_norm = df.copy()
    df_norm['norm_factor'] = df_norm['timepoint'].apply(get_factor)

    for rep in REPS:
        df_norm[rep] = df_norm[rep] / df_norm['norm_factor']

    df_norm.drop(columns=['norm_factor'], inplace=True)

    # Save this initial normalized DataFrame in memory
    normalized_dfs[norm_label] = df_norm

    # Write out the initial normalized table to TSV
    out_path = os.path.join(OUTPUT_DIR, f'{norm_label}_normalized.tsv')
    df_norm.to_csv(out_path, sep='\t', index=False)
    print(f"Wrote normalized file: {out_path}")

# Step 2: Load the 10?~@~Psampling summary and apply each to both normalized_dfs
sampling_df = pd.read_csv(SAMPLING_FILE)
required_sampling_cols = {'Sampling'} | {f'tp{i}' for i in range(1, 17)}
if not required_sampling_cols.issubset(sampling_df.columns):
    missing = required_sampling_cols - set(sampling_df.columns)
    raise ValueError(f"Missing required columns in {SAMPLING_FILE}: {missing}")
  # For each sampling row (Sampling_1 ?~@? Sampling_10):
for _, row in sampling_df.iterrows():
    sampling_name = row['Sampling']  # e.g. "Sampling_1"
    # Extract the index number (1..10) from the name
    try:
        samp_index = int(sampling_name.split('_')[-1])
    except ValueError:
        raise ValueError(f"Unexpected Sampling format: {sampling_name}")

    # Build sampling-specific factors: tpX / tp1 for this row
    samp_factors = {}
    tp1_val = float(row['tp1'])
    if tp1_val == 0:
        raise ZeroDivisionError(f"tp1 is zero in sampling row {sampling_name}, cannot divide.")
    for tp_col in [f'tp{i}' for i in range(1, 17)]:
        try:
            numerator = float(row[tp_col])
            samp_factors[tp_col] = numerator / tp1_val
        except ValueError:
            raise ValueError(f"Non numeric value in column '{tp_col}' of {SAMPLING_FILE}, row {sampling_name}")

    # Helper to map a timepoint to this sampling's factor
    def get_sampling_factor(tp):
        tp_str = str(tp).strip()
        if not tp_str.startswith('tp'):
            tp_str = 'tp' + tp_str
        return samp_factors.get(tp_str, 1.0)

    # Now apply to both previously?~@~Pnormalized DataFrames
    for norm_label, base_df in normalized_dfs.items():
        df_samp = base_df.copy()
        df_samp['norm_factor_samp'] = df_samp['timepoint'].apply(get_sampling_factor)

        for rep in REPS:
            df_samp[rep] = df_samp[rep] / df_samp['norm_factor_samp']

        df_samp.drop(columns=['norm_factor_samp'], inplace=True)

        # Determine output filename per user specification
        if norm_label == 'assigned_reads':
            prefix = 'assigned_read_common_genes_normalized'
        else:  # 'read_depth'
            prefix = 'read_depth_common_genes_normalized'

        filename = f"{prefix}_{samp_index}.tsv"
        out_path = os.path.join(OUTPUT_DIR, filename)
        df_samp.to_csv(out_path, sep='\t', index=False)
        print(f"Wrote sampling normalized file: {out_path}")
