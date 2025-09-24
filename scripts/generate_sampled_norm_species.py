import pandas as pd
import os

def add_sum_column():
    input_file = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/common_genes_tp_least10_assigned_read_normalized/common_genes_total_counts_summary_assigned_reads_normalized.csv"
    output_base = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime/common_genes_tp_least10_assigned_read_normalized/tp_log2FC_common_genes_assigned_reads_normalized"

    # Load the CSV file
    df = pd.read_csv(input_file)

    # Sum columns 2 to 17 (indices 1 to 16)
    df["SUM"] = df.iloc[:, 1:17].sum(axis=1)

    # Sort the DataFrame by the SUM column in descending order
    df = df.sort_values(by="SUM", ascending=False)

    # Divide the data into three equal bins
    total_rows = len(df)
    bin_size, remainder = divmod(total_rows, 3)

    num_bin1 = bin_size
    num_bin2 = bin_size
    num_bin3 = bin_size

    if remainder == 1:
        num_bin2 += 1
    elif remainder == 2:
        num_bin2 += 1
        num_bin3 += 1

    bin1 = df.iloc[:num_bin1]
    bin2 = df.iloc[num_bin1:num_bin1 + num_bin2]
    bin3 = df.iloc[num_bin1 + num_bin2:]

    # Compute median SUM values for each bin
    median_bin1 = bin1["SUM"].median()
    median_bin2 = bin2["SUM"].median()
    median_bin3 = bin3["SUM"].median()

    print(f"Rows in Bin 1: {num_bin1}")
    print(f"Rows in Bin 2: {num_bin2}")
    print(f"Rows in Bin 3: {num_bin3}")

    # Scaling function
    def scale_values(df_subset, category_median, reference_median):
        if category_median > 0:
            scaling_factor = reference_median / category_median
            df_subset.iloc[:, 1:17] = df_subset.iloc[:, 1:17] * scaling_factor
            df_subset["SUM"] = df_subset["SUM"] * scaling_factor
        return df_subset

    # Generate 10 output files with correct selection counts and scaling
    sampled_files = []
    for i in range(10):
        sampled_df = pd.concat([
            bin1.sample(n=min(len(bin1), 3), random_state=i, replace=False) if len(bin1) > 0 else pd.DataFrame(),
            scale_values(bin2.sample(n=min(len(bin2), 3), random_state=i, replace=False), median_bin2, median_bin1) if len(bin2) > 0 else pd.DataFrame(),
            scale_values(bin3.sample(n=min(len(bin3), 3), random_state=i, replace=False), median_bin3, median_bin1) if len(bin3) > 0 else pd.DataFrame()
        ])
        output_file_i = f"{output_base}_sampled_{i+1}.csv"
        sampled_df.to_csv(output_file_i, index=False)
        sampled_files.append(output_file_i)
          # Generate summary
    summary_data = []
    for i in range(10):
        file_path = f"{output_base}_sampled_{i+1}.csv"
        if os.path.exists(file_path):
            sampled_df = pd.read_csv(file_path)
            sum_values = sampled_df.iloc[:, list(range(1, 17)) + [17]].sum().tolist()
            summary_data.append([f"Sampling_{i+1}"] + sum_values)

    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_data, columns=["Sampling"] + [f"tp{i}" for i in range(1, 17)] + ["SUM"])
    summary_output_file = f"{output_base}_summary.csv"
    summary_df.to_csv(summary_output_file, index=False)

if __name__ == "__main__":
    add_sum_column()
