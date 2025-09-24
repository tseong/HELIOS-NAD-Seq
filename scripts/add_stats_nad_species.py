#!/usr/bin/env python3
import os
import pandas as pd
import ast

# ---------------------------
# Settings and input paths
# ---------------------------
base_path   = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_Astart"
input_file  = os.path.join(base_path, "common_nad_genes_across_timepoints.csv")
output_file = os.path.join(base_path, "common_nad_genes_across_timepoints_with_stats.csv")

ALL_TPS = [f"tp{i}" for i in range(1, 17)]

# ---------------------------
# Pre-read normalized files into a dictionary
# tp_data["tpX"] = { geneid: (avg, stderr), ... }
# ---------------------------
tp_data = {}
for tp in ALL_TPS:
    norm_file = os.path.join(base_path, tp, "common_21_assigned_read_normalized", f"normalized_nad_readCount_{tp}.csv")
    if os.path.isfile(norm_file):
        try:
            df_norm = pd.read_csv(norm_file, sep=None, engine='python')
            # Assume first col is Geneid; last two are avg and stderr
            gene_col   = df_norm.columns[0]
            avg_col    = df_norm.columns[-2]
            stderr_col = df_norm.columns[-1]
            tp_data[tp] = {
                str(row[gene_col]): (row[avg_col], row[stderr_col])
                for _, row in df_norm.iterrows()
            }
            print(f"Loaded data for {tp} from {norm_file}")
        except Exception as e:
            print(f"Error loading {norm_file}: {e}")
    else:
        print(f"Warning: file not found for {tp}: {norm_file}")

# ---------------------------
# Read the input common file
# ---------------------------
try:
    df_common = pd.read_csv(input_file, sep=None, engine='python')
except Exception as e:
    print(f"Error reading {input_file}: {e}")
    raise SystemExit(1)

# Require Geneid & TimePoints columns (TimePoints can be ignored for coverage)
if "Geneid" not in df_common.columns or "TimePoints" not in df_common.columns:
    print("Input file must contain 'Geneid' and 'TimePoints' columns")
    raise SystemExit(1)
  # ---------------------------
# Build per-gene stats across ALL timepoints
# ---------------------------
def build_all_tp_stats(row):
    gene = str(row["Geneid"])
    # (Optional) keep original list parsed, but we won?~@~Yt rely on it
    try:
        _ = ast.literal_eval(str(row["TimePoints"]))
    except Exception:
        pass

    # For every tp1..tp16, append "tp, avg, stderr", using NA if missing
    entries = []
    for tp in ALL_TPS:
        if tp in tp_data and gene in tp_data[tp]:
            avg, stderr = tp_data[tp][gene]
        else:
            avg, stderr = "NA", "NA"
        entries.append(f"{tp}, {avg}, {stderr}")

    # Return a stringified list for readability/compatibility
    return str(entries)

# Update column
df_common["TimePoints"] = df_common.apply(build_all_tp_stats, axis=1)

# ---------------------------
# Write output
# ---------------------------
df_common.to_csv(output_file, index=False)
print(f"Saved updated file with stats (all tps) to: {output_file}")
