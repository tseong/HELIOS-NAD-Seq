import pandas as pd
import glob
import os
import re

# Base path where the table files are located
base_path = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_Astart"

# Get all timepoint subfolders (tp1 to tp16)
timepoint_folders = [os.path.join(base_path, f"tp{i}") for i in range(1, 17) if os.path.isdir(os.path.join(base_path, f"tp{i}"))]

# Helper function to extract the "barcode" (bc01, bc02, etc.)
def extract_barcode(filename):
    match = re.search(r'(bc\d+)', filename)
    return match.group(1) if match else None

# Process each timepoint folder
for tp_folder in timepoint_folders:
    # Get all paired and unpaired table files in the current timepoint folder
    paired_files = glob.glob(os.path.join(tp_folder, "*_paired*table"))
    unpaired_files = glob.glob(os.path.join(tp_folder, "*_unpaired*table"))

    # Check if files exist
    if not paired_files or not unpaired_files:
        print(f"Skipping {tp_folder}: Missing either paired or unpaired files.")
        continue

    # Group files by barcode
    paired_map = {}
    unpaired_map = {}

    for file in paired_files:
        barcode = extract_barcode(file)
        if barcode:
            paired_map.setdefault(barcode, []).append(file)

    for file in unpaired_files:
        barcode = extract_barcode(file)
        if barcode:
            unpaired_map.setdefault(barcode, []).append(file)

    # Find barcodes present in both paired and unpaired
    common_barcodes = set(paired_map.keys()).intersection(set(unpaired_map.keys()))

    if not common_barcodes:
        print(f"Skipping {tp_folder}: No matching barcodes found across paired and unpaired files.")
        continue
          # Initialize merged DataFrame
    merged_df = None

    # Process each barcode
    for barcode in sorted(common_barcodes):
        all_files = paired_map[barcode] + unpaired_map[barcode]

        # Initialize DataFrame to sum all timepoints for this barcode
        barcode_df = None

        for file in all_files:
            df = pd.read_csv(file, sep="\t", comment='#', header=0, usecols=[0, 6], names=["Geneid", "count"])

            if barcode_df is None:
                barcode_df = df
            else:
                barcode_df = pd.merge(barcode_df, df, on="Geneid", how="outer", suffixes=('', '_new')).fillna(0)
                barcode_df["count"] = barcode_df["count"] + barcode_df["count_new"]
                barcode_df.drop(columns=["count_new"], inplace=True)

        # Rename final count column to the barcode (bc01, bc02, etc.)
        barcode_df.rename(columns={"count": barcode}, inplace=True)

        # Merge into the main merged_df
        if merged_df is None:
            merged_df = barcode_df
        else:
            merged_df = pd.merge(merged_df, barcode_df, on="Geneid", how="outer").fillna(0)

    # Define output file path for the current timepoint folder
    output_filename = os.path.join(tp_folder, "merged_by_barcode_Astart_readCount.csv")

    # Save merged DataFrame to CSV
    merged_df.to_csv(output_filename, sep=",", index=False)

    print(f"Merged table written to {output_filename}")
