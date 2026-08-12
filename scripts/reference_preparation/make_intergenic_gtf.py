```python
#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from itertools import groupby


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
# Read gene features from GTF
# ------------------------------------------------------------

def read_genes(gtf_path):
    """
    Read gene features from a GTF file.

    Returns:
        header_lines : list of GTF header/comment lines
        genes        : list of gene records
    """

    header_lines = []
    genes = []

    with gtf_path.open("r") as fin:

        for line in fin:

            if line.startswith("#"):
                header_lines.append(line)
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 9:
                continue

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
# Generate intergenic features
# ------------------------------------------------------------

def generate_intergenic_features(genes, extension):
    """
    Generate one intergenic feature per gene.

    For '+' strand genes:
        Extend upstream by up to <extension> bp without
        crossing the previous gene.

    For '-' strand genes:
        Extend downstream by up to <extension> bp without
        crossing the next gene.

    Coordinates remain 1-based as required by GTF format.
    """

    genes.sort(key=lambda x: (x[0], x[1]))

    # Store the next gene start for each gene
    next_start = {}

    for chrom, group in groupby(
        genes,
        key=lambda x: x[0]
    ):

        chromosome_genes = list(group)

        for i, gene in enumerate(chromosome_genes):

            (
                _chrom,
                _start,
                _end,
                _source,
                _score,
                _strand,
                _frame,
                attributes,
            ) = gene

            gene_id = parse_attrs(attributes).get(
                "gene_id"
            )

            if gene_id is None:
                continue

            if i < len(chromosome_genes) - 1:
                next_start[gene_id] = (
                    chromosome_genes[i + 1][1]
                )
            else:
                next_start[gene_id] = None


    prev_end = {}

    intergenic_features = []


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

        gene_id = parse_attrs(attributes).get(
            "gene_id",
            "."
        )


        # ----------------------------------------------------
        # Positive strand
        # ----------------------------------------------------

        if strand == "+":

            intergenic_end = end

            previous_end = prev_end.get(
                chrom,
                0
            )

            intergenic_start = max(
                previous_end + 1,
                start - extension,
                1
            )

            # If previous gene overlaps current gene,
            # do not extend upstream
            if previous_end >= start:
                intergenic_start = start


        # ----------------------------------------------------
        # Negative strand
        # ----------------------------------------------------

        elif strand == "-":

            intergenic_start = start

            following_start = next_start.get(
                gene_id
            )

            if following_start is None:

                intergenic_end = end + extension

            else:

                gap = following_start - end - 1

                if gap <= 0:

                    # Overlapping downstream gene
                    intergenic_end = end

                elif gap <= extension:

                    intergenic_end = (
                        following_start - 1
                    )

                else:

                    intergenic_end = (
                        end + extension
                    )


        # ----------------------------------------------------
        # Unexpected strand value
        # ----------------------------------------------------

        else:

            print(
                f"WARNING: Skipping gene with "
                f"unsupported strand '{strand}': "
                f"{gene_id}"
            )

            prev_end[chrom] = max(
                prev_end.get(chrom, 0),
                end
            )

            continue


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


        # Track furthest gene end on this chromosome
        prev_end[chrom] = max(
            prev_end.get(chrom, 0),
            end
        )


    return intergenic_features


# ------------------------------------------------------------
# Write output GTF
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
                "\t".join(feature) + "\n"
            )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate strand-aware extended gene/intergenic "
            "features from a reference GTF annotation for "
            "HELIOS NAD-Seq analysis."
        )
    )

    parser.add_argument(
        "input_gtf",
        help="Input reference GTF file."
    )

    parser.add_argument(
        "output_gtf",
        help="Output intergenic GTF file."
    )

    parser.add_argument(
        "--extension",
        type=int,
        default=100,
        help=(
            "Maximum number of bases to extend from the "
            "gene boundary. Default: 100 bp."
        )
    )

    args = parser.parse_args()


    # --------------------------------------------------------
    # Validate arguments
    # --------------------------------------------------------

    input_gtf = Path(
        args.input_gtf
    ).expanduser().resolve()

    output_gtf = Path(
        args.output_gtf
    ).expanduser().resolve()


    if not input_gtf.is_file():
        raise FileNotFoundError(
            f"Input GTF does not exist: "
            f"{input_gtf}"
        )


    if args.extension < 0:
        raise ValueError(
            "--extension must be >= 0"
        )


    output_gtf.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    print(
        "========================================"
    )

    print(
        "HELIOS NAD-Seq: Reference preparation"
    )

    print(
        "Generate intergenic GTF"
    )

    print(
        "========================================"
    )

    print(
        f"Input GTF:     {input_gtf}"
    )

    print(
        f"Output GTF:    {output_gtf}"
    )

    print(
        f"Extension:     {args.extension} bp"
    )

    print(
        "========================================"
    )


    # --------------------------------------------------------
    # Read input
    # --------------------------------------------------------

    header_lines, genes = read_genes(
        input_gtf
    )


    if not genes:
        raise RuntimeError(
            "No gene features were found "
            "in the input GTF."
        )


    print(
        f"Gene features found: {len(genes)}"
    )


    # --------------------------------------------------------
    # Generate intergenic annotation
    # --------------------------------------------------------

    features = generate_intergenic_features(
        genes,
        args.extension
    )


    # --------------------------------------------------------
    # Write result
    # --------------------------------------------------------

    write_gtf(
        output_gtf,
        header_lines,
        features
    )


    print(
        f"Intergenic features written: "
        f"{len(features)}"
    )

    print(
        "========================================"
    )

    print(
        "Reference preparation completed."
    )

    print(
        "========================================"
    )


if __name__ == "__main__":
    main()
```
