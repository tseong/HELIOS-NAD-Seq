#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


# ============================================================
# HELIOS NAD-Seq reference preparation
#
# Generate a strand-aware variable 3' gene annotation.
#
# Expected workflow structure:
#
# tp_workflow/
# ├── scripts/
# │   └── reference_prep/
# │       └── make_variable_3prime_gtf.py
# │
# └── gtf/
#     ├── ncbi_dataset/
#     │   └── data/
#     │       └── GCF_000005845.2/
#     │           └── genomic.gtf
#     │
#     └── GCF_000005845.2_genomic_variable_3prime.gtf
#
#
# Default usage:
#
#   python scripts/reference_prep/make_variable_3prime_gtf.py
#
#
# Optional:
#
#   python scripts/reference_prep/make_variable_3prime_gtf.py \
#       --input-gtf /path/to/input.gtf \
#       --output-gtf /path/to/output.gtf
#
# ============================================================


# ------------------------------------------------------------
# GTF attribute parser
# ------------------------------------------------------------

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attrs(attribute_string):
    """
    Parse GTF attributes into a dictionary.
    """

    return {
        key: value
        for key, value in ATTR_RE.findall(attribute_string)
    }


# ------------------------------------------------------------
# Determine workflow root
# ------------------------------------------------------------

def find_workflow_root(script_path):
    """
    Search upward from this script until a directory containing
    'gtf/' is found.

    This allows the script to live at:

        tp_workflow/scripts/reference_prep/

    while correctly identifying:

        tp_workflow/

    as the workflow root.
    """

    current = script_path.resolve().parent

    for candidate in [current, *current.parents]:

        if (candidate / "gtf").is_dir():
            return candidate

    raise RuntimeError(
        "Could not determine the workflow root.\n"
        "No parent directory containing a 'gtf/' directory "
        "was found."
    )


# ------------------------------------------------------------
# Determine whether a gene is tRNA/rRNA
# ------------------------------------------------------------

def is_trna_or_rrna(attributes):
    """
    Determine whether a gene is explicitly annotated as
    tRNA or rRNA.

    Only biotype/type fields are considered.

    This avoids classifying protein-coding genes such as
    'isoleucine--tRNA ligase' as tRNA simply because the word
    'tRNA' appears in a product description.
    """

    attrs = parse_attrs(attributes)

    fields = [
        attrs.get("gene_biotype", ""),
        attrs.get("transcript_biotype", ""),
        attrs.get("gene_type", ""),
        attrs.get("transcript_type", ""),
    ]

    for value in fields:

        normalized = value.strip().lower()

        if normalized in {
            "trna",
            "rrna"
        }:
            return True

    return False


# ============================================================
# Main
# ============================================================

def main():

    # --------------------------------------------------------
    # Determine workflow root
    # --------------------------------------------------------

    script_path = Path(__file__).resolve()

    workflow_dir = find_workflow_root(
        script_path
    )


    # --------------------------------------------------------
    # Default GCF files for this workflow
    # --------------------------------------------------------

    default_input_gtf = (
        workflow_dir
        / "gtf"
        / "ncbi_dataset"
        / "data"
        / "GCF_000005845.2"
        / "genomic.gtf"
    )


    default_output_gtf = (
        workflow_dir
        / "gtf"
        / "GCF_000005845.2_genomic_variable_3prime.gtf"
    )


    # --------------------------------------------------------
    # Command-line arguments
    # --------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=(
            "Generate a strand-aware variable 3' gene annotation "
            "from a GCF_000005845.2 reference GTF for "
            "HELIOS NAD-Seq analysis."
        )
    )


    parser.add_argument(
        "--input-gtf",
        default=str(default_input_gtf),
        help=(
            "Input reference GTF. Default: "
            "<workflow>/gtf/ncbi_dataset/data/"
            "GCF_000005845.2/genomic.gtf"
        )
    )


    parser.add_argument(
        "--output-gtf",
        default=str(default_output_gtf),
        help=(
            "Output variable 3' GTF. Default: "
            "<workflow>/gtf/"
            "GCF_000005845.2_genomic_variable_3prime.gtf"
        )
    )


    parser.add_argument(
        "--short-gene-threshold",
        type=int,
        default=100,
        help=(
            "Genes shorter than this threshold are retained "
            "in full. Default: 100 bp."
        )
    )


    parser.add_argument(
        "--medium-gene-threshold",
        type=int,
        default=200,
        help=(
            "Upper length threshold for genes using the medium "
            "5' offset. Default: 200 bp."
        )
    )


    parser.add_argument(
        "--medium-offset",
        type=int,
        default=50,
        help=(
            "Bases removed from the biological 5' side of genes "
            "between the short and medium thresholds. "
            "Default: 50 bp."
        )
    )


    parser.add_argument(
        "--long-offset",
        type=int,
        default=100,
        help=(
            "Bases removed from the biological 5' side of genes "
            "longer than the medium threshold. "
            "Default: 100 bp."
        )
    )


    args = parser.parse_args()


    # --------------------------------------------------------
    # Resolve paths
    # --------------------------------------------------------

    input_gtf = Path(
        args.input_gtf
    ).expanduser().resolve()


    output_gtf = Path(
        args.output_gtf
    ).expanduser().resolve()


    # --------------------------------------------------------
    # Validate arguments
    # --------------------------------------------------------

    if not input_gtf.is_file():

        raise FileNotFoundError(
            "Input GTF does not exist:\n"
            f"{input_gtf}\n\n"
            "Expected default location:\n"
            f"{default_input_gtf}\n\n"
            "Alternatively specify it explicitly with:\n"
            "--input-gtf /path/to/genomic.gtf"
        )


    if args.short_gene_threshold < 1:

        raise ValueError(
            "--short-gene-threshold must be >= 1."
        )


    if (
        args.medium_gene_threshold
        < args.short_gene_threshold
    ):

        raise ValueError(
            "--medium-gene-threshold must be greater than or "
            "equal to --short-gene-threshold."
        )


    if args.medium_offset < 0:

        raise ValueError(
            "--medium-offset must be >= 0."
        )


    if args.long_offset < 0:

        raise ValueError(
            "--long-offset must be >= 0."
        )


    output_gtf.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Print configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq reference preparation")
    print("Generate variable 3' GTF")
    print("========================================")

    print(
        f"Workflow directory:      {workflow_dir}"
    )

    print(
        f"Input GTF:               {input_gtf}"
    )

    print(
        f"Output GTF:              {output_gtf}"
    )

    print(
        f"Short gene threshold:    "
        f"{args.short_gene_threshold} bp"
    )

    print(
        f"Medium gene threshold:   "
        f"{args.medium_gene_threshold} bp"
    )

    print(
        f"Medium offset:            "
        f"{args.medium_offset} bp"
    )

    print(
        f"Long offset:              "
        f"{args.long_offset} bp"
    )

    print("========================================")


    # --------------------------------------------------------
    # Counters
    # --------------------------------------------------------

    genes_processed = 0

    trna_rrna_genes = 0

    short_genes = 0
    medium_genes = 0
    long_genes = 0

    unsupported_strand = 0
    malformed_lines = 0
    invalid_coordinates = 0


    # --------------------------------------------------------
    # Process GTF
    # --------------------------------------------------------

    with input_gtf.open("r") as gtf, \
         output_gtf.open("w") as out_gtf:

        for line in gtf:

            # Preserve GTF headers/comments
            if line.startswith("#"):

                out_gtf.write(line)

                continue


            columns = line.rstrip("\n").split("\t")


            # Skip malformed records
            if len(columns) < 9:

                malformed_lines += 1

                continue


            # Keep only gene features
            feature_type = columns[2]


            if feature_type != "gene":

                continue


            # ------------------------------------------------
            # Parse coordinates
            # ------------------------------------------------

            try:

                start = int(columns[3])
                end = int(columns[4])

            except ValueError:

                print(
                    "WARNING: Non-integer coordinates detected: "
                    f"{columns[0]}:{columns[3]}-{columns[4]}"
                )

                invalid_coordinates += 1

                continue


            strand = columns[6]

            attributes = columns[8]


            # ------------------------------------------------
            # Validate strand
            # ------------------------------------------------

            if strand not in {"+", "-"}:

                print(
                    f"WARNING: Skipping feature with unsupported "
                    f"strand '{strand}' at "
                    f"{columns[0]}:{start}-{end}"
                )

                unsupported_strand += 1

                continue


            # ------------------------------------------------
            # Validate coordinates
            # ------------------------------------------------

            gene_length = (
                end - start + 1
            )


            if gene_length <= 0:

                print(
                    f"WARNING: Invalid gene coordinates: "
                    f"{columns[0]}:{start}-{end}"
                )

                invalid_coordinates += 1

                continue


            genes_processed += 1


            # ------------------------------------------------
            # Determine whether this is explicitly annotated
            # as tRNA/rRNA
            # ------------------------------------------------

            trna_rrna = is_trna_or_rrna(
                attributes
            )


            # ------------------------------------------------
            # Determine biological 5' offset
            #
            # tRNA/rRNA:
            #     keep complete gene
            #
            # <100 bp:
            #     keep complete gene
            #
            # 100-200 bp:
            #     remove 50 bp from 5' side
            #
            # >200 bp:
            #     remove 100 bp from 5' side
            # ------------------------------------------------

            if trna_rrna:

                offset = 0

                trna_rrna_genes += 1


            elif (
                gene_length
                < args.short_gene_threshold
            ):

                offset = 0

                short_genes += 1


            elif (
                gene_length
                <= args.medium_gene_threshold
            ):

                offset = args.medium_offset

                medium_genes += 1


            else:

                offset = args.long_offset

                long_genes += 1


            # Prevent offset from removing the whole gene
            offset = min(
                offset,
                gene_length - 1
            )


            # ------------------------------------------------
            # Strand-aware removal from biological 5' end
            # ------------------------------------------------

            if strand == "+":

                new_start = (
                    start + offset
                )

                new_end = end


            else:

                new_start = start

                new_end = (
                    end - offset
                )


            # ------------------------------------------------
            # Final coordinate safety check
            # ------------------------------------------------

            if new_start > new_end:

                print(
                    "WARNING: Generated invalid coordinates for "
                    f"{columns[0]}:{start}-{end}. "
                    "Keeping original gene coordinates."
                )

                new_start = start
                new_end = end


            # ------------------------------------------------
            # Update GTF coordinates
            # ------------------------------------------------

            columns[3] = str(
                new_start
            )

            columns[4] = str(
                new_end
            )


            # ------------------------------------------------
            # Write feature
            # ------------------------------------------------

            out_gtf.write(
                "\t".join(columns)
                + "\n"
            )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()

    print("========================================")
    print("Reference preparation completed.")
    print("----------------------------------------")

    print(
        f"Genes processed:          {genes_processed}"
    )

    print(
        f"tRNA/rRNA genes:          {trna_rrna_genes}"
    )

    print(
        f"Short genes:              {short_genes}"
    )

    print(
        f"Medium genes:             {medium_genes}"
    )

    print(
        f"Long genes:               {long_genes}"
    )

    print(
        f"Unsupported strand:       {unsupported_strand}"
    )

    print(
        f"Malformed GTF lines:      {malformed_lines}"
    )

    print(
        f"Invalid coordinates:      {invalid_coordinates}"
    )

    print(
        f"Total output genes:       {genes_processed}"
    )

    print()

    print(
        "Output GTF:"
    )

    print(
        output_gtf
    )

    print("========================================")


if __name__ == "__main__":
    main()
