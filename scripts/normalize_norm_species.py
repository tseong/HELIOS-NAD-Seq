import pandas as pd

# Input files
counts_file = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/common_genes_tp_least10_assigned_read_normalized/common_genes_total_counts_summary.csv"
read_depth_file = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/assigned_reads.csv"
output_file = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/common_genes_tp_least10_assigned_read_normalized/common_genes_total_counts_summary_assigned_reads_normalized.csv"

# Load the counts and read depth data
counts_df = pd.read_csv(counts_file)
read_depth_df = pd.read_csv(read_depth_file)

# Extract normalization factors: tp1 / tpX
tp_columns = [f"tp{i}" for i in range(1, 17)]
read_depths = read_depth_df.loc[0, tp_columns]

norm_factors = read_depths["tp1"] / read_depths

# Normalize the count data
normalized_df = counts_df.copy()
normalized_df[tp_columns] = normalized_df[tp_columns].mul(norm_factors.values, axis=1)

# Calculate %CV across tp1?~@~Stp16 for each gene
normalized_df["%CV"] = (normalized_df[tp_columns].std(axis=1) / normalized_df[tp_columns].mean(axis=1)) * 100

# Save the output
normalized_df.to_csv(output_file, index=False)
print(f"Normalized gene counts with %CV saved to {output_file}")
