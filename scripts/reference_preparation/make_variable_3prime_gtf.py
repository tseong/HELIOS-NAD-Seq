#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate a strand-aware variable 3' gene annotation "
            "from a reference GTF for HELIOS NAD-Seq analysis."
        )
    )

    parser.add_argument(
        "input_gtf",
        help="Input reference GTF file."
    )

    parser.add_argument(
        "output_gtf",
        help="Output variable 3' GTF file."
    )

    parser.add_argument(
        "--short-gene-threshold",
        type=int,
        default=100,
        help=(
            "Genes shorter than this threshold are retained in full. "
            "Default: 100 bp."
        )
    )

    parser.add_argument(
        "--medium-gene-threshold",
        type=int,
        default=200,
        help=(
            "Upper length threshold for genes using the medium "
            "offset. Default: 200 bp."
        )
    )

    parser.add_argument(
        "--medium-offset",
        type=int,
        default=50,
        help=(
            "Number of bases removed from the 5' side of genes "
            "between the short and medium thresholds. "
            "Default: 50 bp."
        )
    )

    parser.add_argument(
        "--long-offset",
        type=int,
        default=100,
        help=(
            "Number of bases removed from the 5' side of genes "
            "longer than the medium threshold. "
            "Default: 100 bp."
        )
    )

    args = parser.parse_args()


    # --------------------------------------------------------
    # Input/output paths
    # --------------------------------------------------------

    input_gtf = Path(
        args.input_gtf
    ).expanduser().resolve()

    output_gtf = Path(
        args.output_gtf
    ).expanduser().resolve()


    if not input_gtf.is_file():
        raise FileNotFoundError(
            f"Input GTF does not exist: {input_gtf}"
        )


    if args.short_gene_threshold < 1:
        raise ValueError(
            "--short-gene-threshold must be >= 1"
        )

    if (
        args.medium_gene_threshold
        < args.short_gene_threshold
    ):
        raise ValueError(
            "--medium-gene-threshold must be greater than or "
            "equal to --short-gene-threshold"
        )

    if args.medium_offset < 0:
        raise ValueError(
            "--medium-offset must be >= 0"
        )

    if args.long_offset < 0:
        raise ValueError(
            "--long-offset must be >= 0"
        )


    output_gtf.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    print("========================================")
    print("HELIOS NAD-Seq: Reference preparation")
    print("Generate variable 3' GTF")
    print("========================================")
    print(f"Input GTF:              {input_gtf}")
    print(f"Output GTF:             {output_gtf}")
    print(
        f"Short gene threshold:   "
        f"{args.short_gene_threshold} bp"
    )
    print(
        f"Medium gene threshold:  "
        f"{args.medium_gene_threshold} bp"
    )
    print(
        f"Medium offset:           "
        f"{args.medium_offset} bp"
    )
    print(
        f"Long offset:             "
        f"{args.long_offset} bp"
    )
    print("========================================")


    # --------------------------------------------------------
    # Process GTF
    # --------------------------------------------------------

    output_lines = []

    genes_processed = 0
    short_genes = 0
    medium_genes = 0
    long_genes = 0
    trna_rrna_genes = 0


    with input_gtf.open("r") as gtf:

        for line in gtf:

            # Preserve original GTF headers/comments
            if line.startswith("#"):
                output_lines.append(line)
                continue


            columns = line.rstrip("\n").split("\t")


            if len(columns) < 9:
                continue


            feature_type = columns[2]


            # Only retain gene features
            if feature_type != "gene":
                continue


            start = int(columns[3])
            end = int(columns[4])
            strand = columns[6]
            attributes = columns[8]


            if strand not in {"+", "-"}:
                print(
                    f"WARNING: Skipping feature with unsupported "
                    f"strand '{strand}' at "
                    f"{columns[0]}:{start}-{end}"
                )
                continue


            gene_length = end - start + 1


            if gene_length <= 0:
                print(
                    f"WARNING: Skipping invalid gene coordinates: "
                    f"{columns[0]}:{start}-{end}"
                )
                continue


            genes_processed += 1


            # ------------------------------------------------
            # Determine whether gene is tRNA/rRNA
            # ------------------------------------------------

            is_trna_rrna = bool(
                re.search(
                    r'transcript_biotype "(tRNA|rRNA)"',
                    attributes
                )
            )


            # ------------------------------------------------
            # Determine 5' offset
            #
            # tRNA/rRNA:
            #     retain full gene
            #
            # gene length < short threshold:
            #     retain full gene
            #
            # short threshold <= length <= medium threshold:
            #     remove medium_offset bp from 5' side
            #
            # length > medium threshold:
            #     remove long_offset bp from 5' side
            # ------------------------------------------------

            if is_trna_rrna:

                offset = 0
                trna_rrna_genes += 1

            elif gene_length < args.short_gene_threshold:

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


            # Prevent an offset from exceeding gene length
            offset = min(
                offset,
                gene_length - 1
            )


            # ------------------------------------------------
            # Apply strand-aware 5' trimming
            # ------------------------------------------------

            if strand == "+":

                new_start = start + offset
                new_end = end

            else:

                new_start = start
                new_end = end - offset


            # ------------------------------------------------
            # Write modified gene coordinates
            # ------------------------------------------------

            columns[3] = str(new_start)
            columns[4] = str(new_end)

            output_lines.append(
                "\t".join(columns) + "\n"
            )


    # --------------------------------------------------------
    # Write output
    # --------------------------------------------------------

    with output_gtf.open("w") as out_gtf:
        out_gtf.writelines(output_lines)


    print()
    print("========================================")
    print("Reference preparation completed.")
    print("----------------------------------------")
    print(f"Genes processed:        {genes_processed}")
    print(f"tRNA/rRNA genes:        {trna_rrna_genes}")
    print(f"Short genes:            {short_genes}")
    print(f"Medium genes:           {medium_genes}")
    print(f"Long genes:             {long_genes}")
    print(f"Output GTF:             {output_gtf}")
    print("========================================")


if __name__ == "__main__":
    main()
