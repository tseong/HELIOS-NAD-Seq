#!/usr/bin/env python3

import argparse
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict


# ============================================================
# HELIOS NAD-Seq pipeline
# SI figure/table preparation
#
# Generate count-weighted sequence logos from FASTA files.
#
# Expected FASTA header:
#
#   >COUNT::CHROM:START-END(STRAND)
#
# Example:
#
#   >25::NC_000913.3:1000-1044(+)
#
# The sequence itself is expected to cover:
#
#   -40 ... TSS(+1) ... +5
#
# i.e. 45 nt total.
#
# Usage:
#
# python scripts/SI_figure_table_prep/make_sequence_logos.py \
#     <WORKFLOW_DIR> \
#     --input-dir PATH
#
# Example:
#
# python scripts/SI_figure_table_prep/make_sequence_logos.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow \
#     --input-dir sequence_logo_input
#
# Default output:
#
#   <WORKFLOW_DIR>/SI_figure_table_prep/sequence_logos/
#
# ============================================================


# ============================================================
# Defaults
# ============================================================

WINDOW_START = -40
WINDOW_END = 4

# -40 through +4 relative to the original zero-based position
# gives 45 bases.
EXPECTED_LEN = (
    WINDOW_END
    - WINDOW_START
    + 1
)

LOGO_FORMAT = "pdf"
DPI = 600


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Generate count-weighted WebLogo sequence logos "
            "for HELIOS NAD-Seq sequence windows."
        )
    )


    parser.add_argument(
        "workflow_dir",
        help="Root tp_workflow directory."
    )


    parser.add_argument(
        "--input-dir",
        required=True,
        help=(
            "Directory containing input FASTA files. "
            "Relative paths are interpreted relative to "
            "WORKFLOW_DIR."
        )
    )


    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. Default: "
            "<WORKFLOW_DIR>/SI_figure_table_prep/"
            "sequence_logos"
        )
    )


    parser.add_argument(
        "--format",
        default="pdf",
        choices=[
            "pdf",
            "png",
            "svg",
            "eps"
        ],
        help="WebLogo output format. Default: pdf."
    )


    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Output resolution. Default: 600."
    )


    parser.add_argument(
        "--window-start",
        type=int,
        default=-40,
        help=(
            "First relative sequence position. "
            "Default: -40."
        )
    )


    parser.add_argument(
        "--window-end",
        type=int,
        default=4,
        help=(
            "Last original relative sequence position. "
            "Default: +4."
        )
    )


    parser.add_argument(
        "--keep-weighted-fasta",
        action="store_true",
        help=(
            "Keep temporary count-expanded FASTA files."
        )
    )


    return parser.parse_args()


# ============================================================
# Generate WebLogo position labels
# ============================================================

def build_position_labels(
    window_start,
    window_end
):
    """
    Convert the original coordinate system:

        -40 ... -1, 0, +1 ... +4

    into TSS notation:

        -40 ... -1, +1, +2 ... +5

    Thus original position 0 becomes TSS +1.
    """

    labels = []


    for pos in range(
        window_start,
        window_end + 1
    ):

        if pos < 0:

            labels.append(
                str(pos)
            )

        else:

            # 0 -> +1
            # 1 -> +2
            # ...
            labels.append(
                f"+{pos + 1}"
            )


    return labels


# ============================================================
# Parse FASTA
# ============================================================

def read_weighted_sequences(
    fasta_path
):
    """
    Read FASTA records and group sequences by chromosome.

    Expected header:

        >COUNT::CHROM:START-END(STRAND)

    COUNT determines the weight of each sequence.
    """

    chrom_seqs = defaultdict(
        list
    )


    current_count = 1
    current_chrom = "unknown"


    with fasta_path.open("r") as fin:

        sequence_parts = []


        def store_sequence():

            if not sequence_parts:
                return


            seq = "".join(
                sequence_parts
            ).strip().upper()


            if not seq:
                return


            chrom_seqs[
                current_chrom
            ].extend(
                [seq]
                * current_count
            )


        for line in fin:

            line = line.strip()


            if not line:
                continue


            # ------------------------------------------------
            # New FASTA header
            # ------------------------------------------------

            if line.startswith(">"):

                # Save previous sequence
                store_sequence()

                sequence_parts = []


                header = line[
                    1:
                ].strip()


                parts = header.split(
                    "::",
                    1
                )


                # --------------------------------------------
                # Sequence weight
                # --------------------------------------------

                try:

                    current_count = int(
                        parts[0]
                    )

                except (
                    ValueError,
                    TypeError
                ):

                    current_count = 1


                if current_count < 1:

                    current_count = 1


                # --------------------------------------------
                # Chromosome
                # --------------------------------------------

                if len(parts) > 1:

                    rest = parts[1]

                    if ":" in rest:

                        current_chrom = (
                            rest.split(
                                ":",
                                1
                            )[0]
                        )

                    else:

                        current_chrom = (
                            rest
                            if rest
                            else "unknown"
                        )

                else:

                    current_chrom = (
                        "unknown"
                    )


            else:

                sequence_parts.append(
                    line
                )


        # Save final FASTA sequence
        store_sequence()


    return chrom_seqs


# ============================================================
# Write weighted FASTA
# ============================================================

def write_weighted_fasta(
    sample,
    chrom,
    seqs,
    output_path
):

    with output_path.open(
        "w"
    ) as fout:

        for i, seq in enumerate(
            seqs,
            start=1
        ):

            fout.write(
                f">{sample}_{chrom}_{i}\n"
            )

            fout.write(
                f"{seq}\n"
            )


# ============================================================
# WebLogo
# ============================================================

def generate_logo(
    weblogo,
    weighted_fasta,
    output_logo,
    sample,
    chrom,
    annotation_string,
    expected_length,
    window_start,
    logo_format,
    dpi
):

    # --------------------------------------------------------
    # Check sequence lengths
    # --------------------------------------------------------

    lengths = set()


    with weighted_fasta.open(
        "r"
    ) as handle:

        for line in handle:

            if line.startswith(">"):
                continue

            seq = line.strip()

            if seq:

                lengths.add(
                    len(seq)
                )


    use_annotation = (
        lengths
        == {
            expected_length
        }
    )


    if not use_annotation:

        print(
            f"  WARNING: {sample} | {chrom}: "
            f"sequence lengths "
            f"{sorted(lengths)} "
            f"!= expected {expected_length}."
        )

        print(
            "  Using natural WebLogo indexing."
        )


    # --------------------------------------------------------
    # WebLogo command
    # --------------------------------------------------------

    cmd = [
        weblogo,

        "-f",
        str(
            weighted_fasta
        ),

        "-o",
        str(
            output_logo
        ),

        "--format",
        logo_format,

        "--resolution",
        str(
            dpi
        ),

        "--title",
        (
            f"{sample} | {chrom}"
        ),

        "--xlabel",
        "Position relative to TSS",

        "--units",
        "bits",

        "--size",
        "large",

        "--number-interval",
        "5",
    ]


    # --------------------------------------------------------
    # Custom TSS labels
    # --------------------------------------------------------

    if use_annotation:

        cmd.extend(
            [
                "--annotate",
                annotation_string
            ]
        )

    else:

        cmd.extend(
            [
                "--first-index",
                str(
                    window_start
                )
            ]
        )


    subprocess.run(
        cmd,
        check=True
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()


    # --------------------------------------------------------
    # Workflow root
    # --------------------------------------------------------

    workflow_dir = Path(
        args.workflow_dir
    ).expanduser().resolve()


    if not workflow_dir.is_dir():

        raise FileNotFoundError(
            f"Workflow directory does not exist:\n"
            f"{workflow_dir}"
        )


    # --------------------------------------------------------
    # Input directory
    # --------------------------------------------------------

    input_dir = Path(
        args.input_dir
    ).expanduser()


    if not input_dir.is_absolute():

        input_dir = (
            workflow_dir
            / input_dir
        )


    input_dir = input_dir.resolve()


    if not input_dir.is_dir():

        raise FileNotFoundError(
            f"Input FASTA directory does not exist:\n"
            f"{input_dir}"
        )


    # --------------------------------------------------------
    # Output directory
    # --------------------------------------------------------

    if args.output_dir:

        output_dir = Path(
            args.output_dir
        ).expanduser()


        if not output_dir.is_absolute():

            output_dir = (
                workflow_dir
                / output_dir
            )


        output_dir = (
            output_dir.resolve()
        )

    else:

        output_dir = (
            workflow_dir
            / "SI_figure_table_prep"
            / "sequence_logos"
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Find WebLogo
    # --------------------------------------------------------

    weblogo = shutil.which(
        "weblogo"
    )


    if weblogo is None:

        raise FileNotFoundError(
            "WebLogo executable was not found in PATH.\n"
            "Activate the Conda environment containing WebLogo "
            "before running this script."
        )


    # --------------------------------------------------------
    # Position information
    # --------------------------------------------------------

    expected_length = (
        args.window_end
        - args.window_start
        + 1
    )


    labels = build_position_labels(
        args.window_start,
        args.window_end
    )


    annotation_string = (
        ",".join(
            labels
        )
    )


    if len(labels) != expected_length:

        raise RuntimeError(
            "Internal position-label length mismatch."
        )


    # --------------------------------------------------------
    # Find input FASTAs
    # --------------------------------------------------------

    fasta_files = sorted(
        list(
            input_dir.glob(
                "*.fa"
            )
        )
        +
        list(
            input_dir.glob(
                "*.fasta"
            )
        )
    )


    if not fasta_files:

        raise FileNotFoundError(
            "No .fa or .fasta files found in:\n"
            f"{input_dir}"
        )


    # --------------------------------------------------------
    # Configuration
    # --------------------------------------------------------

    print("========================================")
    print("HELIOS NAD-Seq")
    print("SI sequence-logo preparation")
    print("========================================")

    print(
        f"Input directory:\n"
        f"  {input_dir}"
    )

    print(
        f"Output directory:\n"
        f"  {output_dir}"
    )

    print(
        f"WebLogo:\n"
        f"  {weblogo}"
    )

    print(
        f"Window: "
        f"{args.window_start} to "
        f"+{args.window_end}"
    )

    print(
        f"Expected sequence length: "
        f"{expected_length}"
    )

    print(
        f"FASTA files found: "
        f"{len(fasta_files)}"
    )

    print("========================================")


    # ========================================================
    # Process FASTAs
    # ========================================================

    total_logos = 0
    failed_logos = 0


    for fasta_path in fasta_files:

        sample = (
            fasta_path.stem
        )


        print()
        print("----------------------------------------")

        print(
            f"Processing: "
            f"{fasta_path.name}"
        )


        # ----------------------------------------------------
        # Parse weighted sequences by chromosome
        # ----------------------------------------------------

        chrom_seqs = (
            read_weighted_sequences(
                fasta_path
            )
        )


        if not chrom_seqs:

            print(
                "  WARNING: No sequences found."
            )

            continue


        # ----------------------------------------------------
        # Process each reference sequence/chromosome
        # ----------------------------------------------------

        for chrom, seqs in (
            sorted(
                chrom_seqs.items()
            )
        ):

            if not seqs:

                continue


            print(
                f"  {chrom}: "
                f"{len(seqs)} weighted sequences"
            )


            # ------------------------------------------------
            # Temporary weighted FASTA
            # ------------------------------------------------

            weighted_fasta = (
                output_dir
                / (
                    f"{sample}_"
                    f"{chrom}_weighted.fa"
                )
            )


            write_weighted_fasta(
                sample,
                chrom,
                seqs,
                weighted_fasta
            )


            # ------------------------------------------------
            # Logo output
            # ------------------------------------------------

            output_logo = (
                output_dir
                / (
                    f"{sample}_"
                    f"{chrom}."
                    f"{args.format}"
                )
            )


            try:

                generate_logo(
                    weblogo=weblogo,
                    weighted_fasta=weighted_fasta,
                    output_logo=output_logo,
                    sample=sample,
                    chrom=chrom,
                    annotation_string=annotation_string,
                    expected_length=expected_length,
                    window_start=args.window_start,
                    logo_format=args.format,
                    dpi=args.dpi
                )


                total_logos += 1


                print(
                    f"  Created: "
                    f"{output_logo.name}"
                )


            except subprocess.CalledProcessError as exc:

                failed_logos += 1


                print(
                    f"  ERROR generating logo "
                    f"for {chrom}: {exc}"
                )


            # ------------------------------------------------
            # Remove temporary FASTA
            # ------------------------------------------------

            if (
                not args.keep_weighted_fasta
                and weighted_fasta.exists()
            ):

                weighted_fasta.unlink()


    # ========================================================
    # Finish
    # ========================================================

    print()
    print("========================================")
    print("Sequence-logo preparation completed")
    print("----------------------------------------")

    print(
        f"Logos generated: "
        f"{total_logos}"
    )

    print(
        f"Failed:          "
        f"{failed_logos}"
    )

    print()

    print(
        f"Output directory:\n"
        f"{output_dir}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
