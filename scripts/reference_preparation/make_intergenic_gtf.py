#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from itertools import groupby


# ============================================================
# HELIOS NAD-Seq reference preparation
#
# Generate strand-aware extended gene features for
# featureCounts.
#
# Expected workflow structure:
#
# tp_workflow/
# ├── scripts/
# │   └── reference_prep/
# │       └── make_intergenic_gtf.py
# │
# └── gtf/
#     ├── ncbi_dataset/
#     │   └── data/
#     │       └── GCF_000005845.2/
#     │           └── genomic.gtf
#     │
#     └── GCF_000005845.2_genomic_intergenic.gtf
#
#
# Usage:
#
#   python scripts/reference_prep/make_intergenic_gtf.py
#
#
# Optional:
#
#   python scripts/reference_prep/make_intergenic_gtf.py \
#       --input-gtf /path/to/input.gtf \
#       --output-gtf /path/to/output.gtf \
#       --extension 100
#
# ============================================================


# ------------------------------------------------------------
# GTF attribute parsing
# ------------------------------------------------------------

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')


def parse_attrs(attribute_string):
    """
    Parse a GTF attribute string into a dictionary.
    """

    return {
        key: value
        for key, value in ATTR_RE.findall(attribute_string)
    }


# ------------------------------------------------------------
# Find workflow root
# ------------------------------------------------------------

def find_workflow_root(script_path):
    """
    Search upward from the script location until a directory
    containing 'gtf/' is found.

    This allows the script to live at:

        tp_workflow/scripts/reference_prep/

    while correctly identifying:

        tp_workflow/

    as the workflow root.
    """

    current = script_path.resolve().parent

    for candidate in [current, *current.parents]:

        gtf_dir = candidate / "gtf"

        if gtf_dir.is_dir():
            return candidate

    raise RuntimeError(
        "Could not determine workflow root. "
        "No parent directory containing 'gtf/' was found."
    )


# ------------------------------------------------------------
# Read gene features
# ------------------------------------------------------------

def read_genes(gtf_path):
    """
    Read gene features from a GTF file.

    Returns
    -------
    header_lines : list
        Original GTF comment/header lines.

    genes : list
        Gene records containing chromosome, coordinates,
        source, score, strand, frame and attributes.
    """

    header_lines = []
    genes = []

    with gtf_path.open("r") as fin:

        for line in fin:

            # Preserve GTF headers/comments
            if line.startswith("#"):
                header_lines.append(line)
                continue

            cols = line.rstrip("\n").split("\t")

            # Skip malformed lines
            if len(cols) < 9:
                continue

            # Only gene features are used
            if cols[2] != "gene":
                continue

            chrom = cols[0]
            source = cols[1]
            start = int(cols[3])
            end = int(cols[4])
            score = cols[5]
            strand = cols[6]
            frame = cols[7]
            attributes = cols[8]

            genes.append(
                (
                    chrom,
                    start,
                    end,
                    source,
                    score,
                    strand,
                    frame,
                    attributes,
                )
            )

    return header_lines, genes


# ------------------------------------------------------------
# Generate strand-aware extended features
# ------------------------------------------------------------

def generate_intergenic_features(genes, extension):
    """
    Generate one strand-aware extended feature per gene.

    Positive-strand genes:
        Extend upstream by up to `extension` bp without
        crossing the preceding annotated gene.

    Negative-strand genes:
        Extend downstream by up to `extension` bp without
        crossing the following annotated gene.

    The full annotated gene itself is retained.

    Coordinates remain 1-based and inclusive as required by GTF.
    """

    # Sort by chromosome and genomic position
    genes = sorted(
        genes,
        key=lambda x: (x[0], x[1], x[2])
    )


    # --------------------------------------------------------
    # Determine next gene start for each record
    # --------------------------------------------------------

    next_start = {}

    for chrom, group in groupby(
        genes,
        key=lambda x: x[0]
    ):

        chromosome_genes = list(group)

        for i, gene in enumerate(chromosome_genes):

            record_key = (
                gene[0],
                gene[1],
                gene[2],
                gene[5],
                gene[7]
            )

            if i < len(chromosome_genes) - 1:

                next_start[record_key] = (
                    chromosome_genes[i + 1][1]
                )

            else:

                next_start[record_key] = None


    # --------------------------------------------------------
    # Generate extended features
    # --------------------------------------------------------

    intergenic_features = []

    prev_end = {}

    positive_count = 0
    negative_count = 0
    skipped_count = 0


    for (
        chrom,
        start,
        end,
        source,
        score,
        strand,
        frame,
        attributes,
    ) in genes:

        attrs = parse_attrs(
            attributes
        )

        gene_id = attrs.get(
            "gene_id",
            attrs.get(
                "gene",
                f"{chrom}:{start}-{end}"
            )
        )


        # ----------------------------------------------------
        # Positive strand
        # ----------------------------------------------------

        if strand == "+":

            previous_end = prev_end.get(
                chrom,
                0
            )

            intergenic_start = max(
                1,
                start - extension,
                previous_end + 1
            )

            intergenic_end = end


            # If previous gene overlaps current gene,
            # do not extend upstream
            if previous_end >= start:

                intergenic_start = start


            positive_count += 1


        # ----------------------------------------------------
        # Negative strand
        # ----------------------------------------------------

        elif strand == "-":

            intergenic_start = start

            record_key = (
                chrom,
                start,
                end,
                strand,
                attributes
            )

            following_start = next_start.get(
                record_key
            )


            if following_start is None:

                intergenic_end = (
                    end + extension
                )

            else:

                gap = (
                    following_start
                    - end
                    - 1
                )


                if gap <= 0:

                    # Next gene overlaps current gene
                    intergenic_end = end

                elif gap <= extension:

                    # Extend only until immediately before
                    # the next gene
                    intergenic_end = (
                        following_start - 1
                    )

                else:

                    intergenic_end = (
                        end + extension
                    )


            negative_count += 1


        # ----------------------------------------------------
        # Unsupported strand
        # ----------------------------------------------------

        else:

            print(
                f"WARNING: Skipping gene with unsupported "
                f"strand '{strand}': {gene_id}"
            )

            prev_end[chrom] = max(
                prev_end.get(chrom, 0),
                end
            )

            skipped_count += 1

            continue


        # ----------------------------------------------------
        # Add intergenic feature
        #
        # Step 06 featureCounts uses:
        #
        #   -t intergenic
        #
        # therefore the feature type remains "intergenic".
        # ----------------------------------------------------

        intergenic_features.append(
            [
                chrom,
                source,
                "intergenic",
                str(intergenic_start),
                str(intergenic_end),
                score,
                strand,
                frame,
                attributes,
            ]
        )


        # Track furthest annotated gene end
        prev_end[chrom] = max(
            prev_end.get(chrom, 0),
            end
        )


    return (
        intergenic_features,
        positive_count,
        negative_count,
        skipped_count
    )


# ------------------------------------------------------------
# Write GTF
# ------------------------------------------------------------

def write_gtf(
    output_path,
    header_lines,
    features
):

    with output_path.open("w") as fout:

        for line in header_lines:
            fout.write(line)

        for feature in features:

            fout.write(
                "\t".join(feature)
                + "\n"
            )


# ============================================================
# Main
# ============================================================

def main():

    # --------------------------------------------------------
    # Determine workflow root automatically
    # --------------------------------------------------------

    script_path = Path(
        __file__
    ).resolve()

    workflow_dir = find_workflow_root(
        script_path
    )


    # --------------------------------------------------------
    # Default GCF reference paths
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
        / "GCF_000005845.2_genomic_intergenic.gtf"
    )


    # --------------------------------------------------------
    # Arguments
    # --------------------------------------------------------

    parser = argparse.ArgumentParser(
        description=(
            "Prepare a strand-aware extended E. coli GTF "
            "annotation for HELIOS NAD-Seq featureCounts."
        )
    )


    parser.add_argument(
        "--input-gtf",
        default=str(
            default_input_gtf
        ),
        help=(
            "Input reference GTF. "
            "Default: <workflow>/gtf/ncbi_dataset/data/"
            "GCF_000005845.2/genomic.gtf"
        )
    )


    parser.add_argument(
        "--output-gtf",
        default=str(
            default_output_gtf
        ),
        help=(
            "Output GTF containing extended 'intergenic' "
            "features. Default: "
            "<workflow>/gtf/"
            "GCF_000005845.2_genomic_intergenic.gtf"
        )
    )


    parser.add_argument(
        "--extension",
        type=int,
        default=100,
        help=(
            "Maximum strand-aware extension from the biological "
            "5' gene boundary. Default: 100 bp."
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
    # Validate inputs
    # --------------------------------------------------------

    if not input_gtf.is_file():

        raise FileNotFoundError(
            "Input GTF does not exist:\n"
            f"{input_gtf}\n\n"
            "Expected default location:\n"
            f"{default_input_gtf}\n\n"
            "You can also specify the file manually with:\n"
            "--input-gtf /path/to/genomic.gtf"
        )


    if args.extension < 0:

        raise ValueError(
            "--extension must be >= 0"
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
    print("Generate extended intergenic GTF")
    print("========================================")

    print(
        f"Workflow directory: {workflow_dir}"
    )

    print(
        f"Input GTF:          {input_gtf}"
    )

    print(
        f"Output GTF:         {output_gtf}"
    )

    print(
        f"Extension:          {args.extension} bp"
    )

    print("========================================")


    # --------------------------------------------------------
    # Read GTF
    # --------------------------------------------------------

    header_lines, genes = read_genes(
        input_gtf
    )


    if not genes:

        raise RuntimeError(
            "No 'gene' features were found in the input GTF."
        )


    print(
        f"Gene features found: {len(genes)}"
    )


    # --------------------------------------------------------
    # Generate intergenic features
    # --------------------------------------------------------

    (
        features,
        positive_count,
        negative_count,
        skipped_count
    ) = generate_intergenic_features(
        genes,
        args.extension
    )


    # --------------------------------------------------------
    # Write output GTF
    # --------------------------------------------------------

    write_gtf(
        output_gtf,
        header_lines,
        features
    )


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()

    print("========================================")
    print("Reference preparation completed.")
    print("----------------------------------------")

    print(
        f"Input gene features:     {len(genes)}"
    )

    print(
        f"Positive-strand genes:   {positive_count}"
    )

    print(
        f"Negative-strand genes:   {negative_count}"
    )

    print(
        f"Unsupported/skipped:     {skipped_count}"
    )

    print(
        f"Features written:        {len(features)}"
    )

    print()

    print(
        "Output GTF:"
    )

    print(
        f"{output_gtf}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
