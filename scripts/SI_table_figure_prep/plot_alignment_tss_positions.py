#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from collections import Counter

import matplotlib.pyplot as plt


# ============================================================
# HELIOS NAD-Seq pipeline
# SI figure/table preparation
#
# Plot genomic alignment-start distributions from Step 04 SAMs.
#
# Input:
#
#   <WORKFLOW_DIR>/alignment/
#       <sample>/
#           *_eColi_paired.sam
#           *_eColi_R1_unpaired.sam
#           *_eColi_R2_unpaired.sam
#
# Output:
#
#   <WORKFLOW_DIR>/SI_table_figure_prep/tss_alignment_positions/
#       tp1/
#           bc01/
#           bc02/
#           ...
#
# Usage:
#
# python scripts/SI_table_figure_prep/plot_alignment_tss_positions.py \
#     <WORKFLOW_DIR>
#
# Example:
#
# python scripts/SI_table_figure_prep/plot_alignment_tss_positions.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow
#
# Optional:
#
#   --alignment-dir PATH
#   --output-dir PATH
#   --start-tp 1
#   --end-tp 16
#
# ============================================================


# ============================================================
# Parse SAM and count alignment starts
# ============================================================

def build_tss_counts(sam_files):
    """
    Build per-reference alignment-start counts.

    Returns:
        {
            rname: Counter({
                position: count,
                ...
            })
        }

    SAM POS is 1-based.
    """

    tss = {}


    for sam_path in sam_files:

        with sam_path.open(
            "r"
        ) as handle:

            for line in handle:

                if line.startswith("@"):
                    continue


                fields = line.rstrip(
                    "\n"
                ).split(
                    "\t"
                )


                if len(fields) < 11:
                    continue


                # ------------------------------------------------
                # SAM columns
                # ------------------------------------------------

                try:
                    flag = int(
                        fields[1]
                    )

                    pos = int(
                        fields[3]
                    )

                except ValueError:
                    continue


                rname = fields[2]


                # ------------------------------------------------
                # Skip unmapped reads
                #
                # FLAG 0x4 = unmapped
                # RNAME "*" also indicates unmapped
                # ------------------------------------------------

                if (
                    flag & 0x4
                    or rname == "*"
                    or pos <= 0
                ):

                    continue


                counts = tss.setdefault(
                    rname,
                    Counter()
                )


                counts[
                    pos
                ] += 1


    return tss


# ============================================================
# Plot
# ============================================================

def plot_tss(
    tss_counts,
    output_dir,
    bc,
    tp
):
    """
    Plot alignment-start distribution for each reference.
    """

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    plots_created = 0


    for rname, counts in sorted(
        tss_counts.items()
    ):

        if not counts:
            continue


        total = sum(
            counts.values()
        )


        first = counts.get(
            1,
            0
        )


        # Same criterion as original script:
        # >20% of reads at position 1
        majority_at_first = (
            first > total / 5
        )


        positions = sorted(
            counts
        )


        values = [
            counts[
                pos
            ]
            for pos in positions
        ]


        # ----------------------------------------------------
        # Plot
        # ----------------------------------------------------

        fig, ax = plt.subplots(
            figsize=(8, 3)
        )


        ax.bar(
            positions,
            values
        )


        ax.set_title(
            f"Alignment-start distribution for "
            f"{rname} ({bc}, {tp})"
        )


        ax.set_xlabel(
            "Reference position (1-based)"
        )


        ax.set_ylabel(
            "Read count"
        )


        # ----------------------------------------------------
        # Mark position 1
        # ----------------------------------------------------

        ax.axvline(
            1,
            linestyle="--",
            linewidth=1
        )


        # ----------------------------------------------------
        # Output filename
        # ----------------------------------------------------

        safe_rname = re.sub(
            r"[^A-Za-z0-9_.-]+",
            "_",
            rname
        )


        if majority_at_first:

            filename = (
                f"{safe_rname}_"
                f"tss_majority_position1.pdf"
            )

        else:

            filename = (
                f"{safe_rname}_"
                f"tss_distribution.pdf"
            )


        output_path = (
            output_dir
            / filename
        )


        fig.tight_layout()


        fig.savefig(
            output_path,
            dpi=300,
            bbox_inches="tight"
        )


        plt.close(
            fig
        )


        plots_created += 1


        print(
            f"    {rname}: "
            f"{total} mapped reads; "
            f"position 1 = {first} "
            f"({first / total * 100:.2f}%)"
        )


        print(
            f"      -> {output_path.name}"
        )


    return plots_created


# ============================================================
# Find SAM files for sample
# ============================================================

def find_sample_sams(
    alignment_dir,
    bc,
    tp
):
    """
    Step 04 outputs are stored recursively under sample folders.

    Expected examples:

        bc01_tp1_eColi_paired.sam
        bc01_tp1_eColi_R1_unpaired.sam
        bc01_tp1_eColi_R2_unpaired.sam

    Filename order may contain additional text, so matching
    is done using barcode + timepoint + eColi.
    """

    candidates = []


    for sam_path in alignment_dir.rglob(
        "*_eColi_*.sam"
    ):

        name = sam_path.name


        # Exact barcode token
        bc_match = re.search(
            rf"(^|_){re.escape(bc)}(_|$)",
            name
        )


        # Exact timepoint token, avoiding tp1 matching tp10
        tp_match = re.search(
            rf"(^|_){re.escape(tp)}(_|$)",
            name
        )


        if (
            bc_match
            and tp_match
        ):

            candidates.append(
                sam_path
            )


    return sorted(
        candidates
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Plot alignment-start distributions for "
            "HELIOS NAD-Seq E. coli SAM files."
        )
    )


    parser.add_argument(
        "workflow_dir",
        help="Root tp_workflow directory."
    )


    parser.add_argument(
        "--alignment-dir",
        default=None,
        help=(
            "Step 04 alignment directory. Default: "
            "<WORKFLOW_DIR>/alignment"
        )
    )


    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. Default: "
            "<WORKFLOW_DIR>/SI_figure_table_prep/"
            "tss_alignment_positions"
        )
    )


    parser.add_argument(
        "--start-tp",
        type=int,
        default=1,
        help="First time point. Default: 1."
    )


    parser.add_argument(
        "--end-tp",
        type=int,
        default=16,
        help="Last time point. Default: 16."
    )


    args = parser.parse_args()


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
    # Alignment directory
    # --------------------------------------------------------

    if args.alignment_dir:

        alignment_dir = Path(
            args.alignment_dir
        ).expanduser()


        if not alignment_dir.is_absolute():

            alignment_dir = (
                workflow_dir
                / alignment_dir
            )


        alignment_dir = alignment_dir.resolve()

    else:

        alignment_dir = (
            workflow_dir
            / "alignment"
        )


    if not alignment_dir.is_dir():

        raise FileNotFoundError(
            f"Alignment directory does not exist:\n"
            f"{alignment_dir}"
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


        output_dir = output_dir.resolve()

    else:

        output_dir = (
            workflow_dir
            / "SI_figure_table_prep"
            / "tss_alignment_positions"
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Barcode/timepoint definitions
    # --------------------------------------------------------

    barcodes = [
        f"bc{i:02d}"
        for i in range(
            1,
            9
        )
    ]


    timepoints = [
        f"tp{i}"
        for i in range(
            args.start_tp,
            args.end_tp + 1
        )
    ]


    print("========================================")
    print("HELIOS NAD-Seq")
    print("SI alignment-start distributions")
    print("========================================")

    print(
        f"Alignment directory:\n"
        f"  {alignment_dir}"
    )

    print(
        f"Output directory:\n"
        f"  {output_dir}"
    )

    print(
        f"Time points: "
        f"{timepoints[0]}-{timepoints[-1]}"
    )

    print(
        "Barcodes: bc01-bc08"
    )

    print("========================================")


    total_plots = 0
    samples_processed = 0


    # ========================================================
    # Process each timepoint / barcode
    # ========================================================

    for tp in timepoints:

        for bc in barcodes:

            sam_files = find_sample_sams(
                alignment_dir,
                bc,
                tp
            )


            if not sam_files:
                continue


            print()
            print("----------------------------------------")

            print(
                f"Processing {bc} {tp}: "
                f"{len(sam_files)} SAM files"
            )


            for sam_file in sam_files:

                print(
                    f"  {sam_file}"
                )


            # ------------------------------------------------
            # Combine paired + singleton alignments
            # ------------------------------------------------

            tss_counts = build_tss_counts(
                sam_files
            )


            if not tss_counts:

                print(
                    "  No mapped reads found."
                )

                continue


            sample_output_dir = (
                output_dir
                / tp
                / bc
            )


            created = plot_tss(
                tss_counts,
                sample_output_dir,
                bc,
                tp
            )


            total_plots += created
            samples_processed += 1


    # ========================================================
    # Finish
    # ========================================================

    print()
    print("========================================")
    print("TSS alignment plotting completed")
    print("----------------------------------------")

    print(
        f"Samples processed: "
        f"{samples_processed}"
    )

    print(
        f"Plots generated:   "
        f"{total_plots}"
    )

    print()

    print(
        f"Output directory:\n"
        f"{output_dir}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
