import pandas as pd
from pathlib import Path

# Base paths
base_path       = Path("/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/table_Astart")
nad_genes_file  = base_path / "common_nad_genes_across_timepoints.csv"

# Load NAD gene list (use all Geneid values)
nad_df = pd.read_csv(nad_genes_file, dtype=str)

# Build the global set of genes to extract for every timepoint
all_genes = (
    nad_df["Geneid"]
    .dropna()
    .map(lambda x: x.strip())
    .tolist()
)
all_genes = sorted(set(all_genes))

print(f"Total unique genes from common file: {len(all_genes)}")

# Loop through each timepoint and filter its counts CSV by the global gene list
for i in range(1, 17):
    tp          = f"tp{i}"
    tp_dir      = base_path / tp
    input_file  = tp_dir / "merged_by_barcode_Astart_readCount.csv"
    output_file = tp_dir / "nad_genes_readCount.csv"

    if not input_file.exists():
        print(f"Warning: {input_file} not found. Skipping {tp}.")
        continue

    # Load the counts (must contain a 'Geneid' column)
    df = pd.read_csv(input_file, dtype=str)
    if "Geneid" not in df.columns:
        print(f"Warning: 'Geneid' column not found in {input_file}. Skipping {tp}.")
        continue

    # Normalize Geneid formatting in counts file to match (strip spaces)
    df["Geneid"] = df["Geneid"].astype(str).map(lambda x: x.strip())

    # Keep rows whose Geneid is in the global NAD list
    filtered = df[df["Geneid"].isin(all_genes)].copy()

    # (Optional) sort by Geneid to have consistent ordering
    filtered = filtered.sort_values("Geneid")

    # Write out
    filtered.to_csv(output_file, index=False)
    print(f"{tp}: wrote {len(filtered)} rows to {output_file} (from {len(all_genes)} target genes)")
