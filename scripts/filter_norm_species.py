import os
import pandas as pd

def filter_csv(input_file, output_file):
    """Filter CSV based on log2FoldChange and baseMean values and save to a new file."""
    df = pd.read_csv(input_file)

    # Ensure required columns exist
    if 'log2FoldChange' not in df.columns or 'padj' not in df.columns or 'baseMean' not in df.columns:
        raise ValueError(f"The input CSV {input_file} must contain 'log2FoldChange', 'padj', and 'baseMean' columns.")

    # Apply filtering conditions
    filtered_df = df[(df['log2FoldChange'] >= -1) & (df['log2FoldChange'] <= 1) & (df['baseMean'] >= 100)]

    # Save filtered data to a new CSV file
    filtered_df.to_csv(output_file, index=False)
    print(f"Filtered data saved to {output_file}")

def process_folders(base_path):
    """Process subfolders tp1 to tp16, filtering relevant CSV files."""
    for i in range(1, 17):
        subfolder = f"tp{i}"
        subfolder_path = os.path.join(base_path, subfolder)

        if not os.path.isdir(subfolder_path):
            print(f"Skipping missing folder: {subfolder_path}")
            continue

        # Find the input file that starts with "pydeseq2"
        input_file = None
        for file in os.listdir(subfolder_path):
            if file.startswith("pydeseq2") and file.endswith(".csv"):
                input_file = os.path.join(subfolder_path, file)
                break

        if input_file:
            output_file = os.path.join(subfolder_path, "pydeseq2_results_log2FC_broad_filtered.csv")
            filter_csv(input_file, output_file)
        else:
            print(f"No input file found in {subfolder_path}")

if __name__ == "__main__":
    base_path = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_3prime"
    process_folders(base_path)
