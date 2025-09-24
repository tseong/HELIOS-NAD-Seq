import re

# Define file paths
gtf_file = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/data/GCF_000005845.2_ASM584v2_genomic.gtf"
output_gtf = "/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/data/GCF_000005845.2_ASM584v2_genomic_variable_3prime.gtf"

# Read and process GTF file
filtered_lines = []

with open(gtf_file, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue  # Skip comments

        columns = line.strip().split("\t")
        if len(columns) < 9:
            continue  # Skip malformed lines

        feature_type = columns[2]
        if feature_type != "gene":
            continue  # Keep only "gene" features

        start = int(columns[3])
        end = int(columns[4])
        strand = columns[6]
        attributes = columns[8]
        gene_length = end - start + 1

        # Check if the gene is a tRNA or rRNA
        if re.search(r'transcript_biotype "(tRNA|rRNA)"', attributes):
            window_size = gene_length  # Keep the entire gene
        elif gene_length < 100:
            window_size = gene_length  # Keep the entire gene
        elif 100 <= gene_length <= 200:
            window_size = 50 + (end - start)  # 50 bp downstream of TSS until the end
        else:
            window_size = 100 + (end - start)  # 100 bp downstream of TSS until the end

        if strand == "+":
            if gene_length > 100:  # Only shift for genes > 100 bp
                new_start = start + (50 if gene_length <= 200 else 100)  # Move start downstream
            else:
                new_start = start  # Keep entire gene
            new_end = end
        else:  # Reverse strand (3' is already at end)
            new_start = end - window_size + 1
            new_end = end

        # Modify columns to reflect the new 3' window region
        columns[3] = str(new_start)
        columns[4] = str(new_end)

        filtered_lines.append("\t".join(columns) + "\n")

# Write new GTF file
with open(output_gtf, "w") as out_gtf:
    out_gtf.writelines(filtered_lines)

print(f"Filtered 3' window GTF saved to {output_gtf}")
