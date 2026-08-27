#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# HELIOS NAD-Seq pipeline
# SI figure/table preparation
#
# Plot read-start distributions relative to TSS.
#
# Expected input filenames:
#
#   read_starts_relative_position_bc01_eColi_tp1.bed
#   read_starts_relative_position_bc02_eColi_tp1.bed
#   ...
#
# Expected columns:
#
#   chrom
#   start
#   end
#   rel_pos
#   count
#   strand
#
# Usage:
#
# python scripts/SI_table_figure_prep/plot_tss_read_start_distributions.py \
#     <WORKFLOW_DIR> \
#     --input-dir PATH
#
# Example:
#
# python scripts/SI_table_figure_prep/plot_tss_read_start_distributions.py \
#     /gpfs/bwfor/work/ws/hd_qp355-tae_temp/tp_workflow \
#     --input-dir sam/intermediate
#
# Default output:
#
#   <WORKFLOW_DIR>/SI_figure_table_prep/tss_positions/
#
# ============================================================


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Plot HELIOS read-start distributions relative to TSS "
            "for each barcode, time point and reference sequence."
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
            "Directory containing read_starts_relative_position_*.bed files. "
            "Relative paths are interpreted relative to WORKFLOW_DIR."
        )
    )

    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. Default: "
            "<WORKFLOW_DIR>/SI_figure_table_prep/tss_positions"
        )
    )

    parser.add_argument(
        "--xmin",
        type=int,
        default=-30,
        help="Minimum relative position shown. Default: -30."
    )

    parser.add_argument(
        "--xmax",
        type=int,
        default=30,
        help="Maximum relative position shown. Default: 30."
    )

    parser.add_argument(
        "--tick-step",
        type=int,
        default=5,
        help="X-axis tick interval. Default: 5."
    )

    parser.add_argument(
        "--chroms",
        nargs="*",
        default=None,
        help=(
            "Optional list of reference sequence names to plot. "
            "If omitted, all chromosomes/references found in each BED file "
            "are plotted."
        )
    )

    return parser.parse_args()


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()


    # --------------------------------------------------------
    # Workflow directory
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
            f"Input directory does not exist:\n"
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


        output_dir = output_dir.resolve()

    else:

        output_dir = (
            workflow_dir
            / "SI_figure_table_prep"
            / "tss_positions"
        )


    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Filename pattern
    # --------------------------------------------------------

    pattern = re.compile(
        r"read_starts_relative_position_"
        r"bc(\d{2})_eColi_tp(\d{1,2})\.bed$"
    )


    bed_files = sorted(
        input_dir.glob(
            "read_starts_relative_position_bc*_eColi_tp*.bed"
        )
    )


    if not bed_files:

        raise FileNotFoundError(
            "No read-start BED files found in:\n"
            f"{input_dir}"
        )


    print("========================================")
    print("HELIOS NAD-Seq")
    print("SI TSS read-start plots")
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
        f"BED files found: "
        f"{len(bed_files)}"
    )

    print(
        f"X-axis range: "
        f"{args.xmin} to {args.xmax}"
    )

    if args.chroms:

        print(
            "Reference sequences: "
            + ", ".join(args.chroms)
        )

    else:

        print(
            "Reference sequences: "
            "all references found in each BED file"
        )

    print("========================================")


    plots_created = 0


    # ========================================================
    # Process BED files
    # ========================================================

    for bed_path in bed_files:

        match = pattern.search(
            bed_path.name
        )


        if not match:

            print(
                f"Skipping unrecognized filename: "
                f"{bed_path.name}"
            )

            continue


        bc, tp = match.groups()


        prefix = (
            f"bc{bc}_tp{tp}"
        )


        print()
        print(
            f"Processing: "
            f"{bed_path.name} -> {prefix}"
        )


        # ----------------------------------------------------
        # Load BED-like table
        # ----------------------------------------------------

        df = pd.read_csv(
            bed_path,
            sep="\t",
            header=0
        )


        required_columns = {
            "chrom",
            "start",
            "end",
            "rel_pos",
            "count",
            "strand",
        }


        missing = (
            required_columns
            - set(df.columns)
        )


        if missing:

            raise ValueError(
                f"{bed_path} is missing columns: "
                + ", ".join(
                    sorted(missing)
                )
            )


        df["rel_pos"] = pd.to_numeric(
            df["rel_pos"],
            errors="raise"
        )


        df["count"] = pd.to_numeric(
            df["count"],
            errors="raise"
        )


        # ----------------------------------------------------
        # Determine references to plot
        # ----------------------------------------------------

        if args.chroms:

            chroms = args.chroms

        else:

            chroms = sorted(
                df[
                    "chrom"
                ]
                .dropna()
                .astype(str)
                .unique()
            )


        # ====================================================
        # Plot each reference separately
        # ====================================================

        for chrom in chroms:

            sub = df[
                df[
                    "chrom"
                ].astype(str)
                == str(chrom)
            ].copy()


            if sub.empty:

                continue


            # ------------------------------------------------
            # Aggregate counts by relative position
            #
            # Important if multiple genomic positions share
            # the same relative TSS coordinate.
            # ------------------------------------------------

            sub = (
                sub.groupby(
                    "rel_pos",
                    as_index=False
                )[
                    "count"
                ]
                .sum()
                .sort_values(
                    "rel_pos"
                )
            )


            # ------------------------------------------------
            # Limit plotting window
            # ------------------------------------------------

            window = sub[
                (
                    sub[
                        "rel_pos"
                    ] >= args.xmin
                )
                &
                (
                    sub[
                        "rel_pos"
                    ] <= args.xmax
                )
            ].copy()


            if window.empty:

                print(
                    f"  {chrom}: no reads in "
                    f"{args.xmin}..{args.xmax}; skipping."
                )

                continue


            # ------------------------------------------------
            # Plot
            # ------------------------------------------------

            fig, ax = plt.subplots(
                figsize=(10, 4)
            )


            ax.bar(
                window[
                    "rel_pos"
                ],
                window[
                    "count"
                ],
                width=0.8
            )


            ax.set_xlabel(
                "Position relative to TSS"
            )


            ax.set_ylabel(
                "Read count"
            )


            ax.set_title(
                f"{chrom} read-start distribution: "
                f"{prefix}"
            )


            # ------------------------------------------------
            # TSS marker
            # ------------------------------------------------

            ax.axvline(
                0,
                linewidth=1,
                linestyle="--"
            )


            # ------------------------------------------------
            # X axis
            # ------------------------------------------------

            ticks = list(
                range(
                    args.xmin,
                    args.xmax + 1,
                    args.tick_step
                )
            )


            ax.set_xticks(
                ticks
            )


            ax.set_xlim(
                args.xmin,
                args.xmax
            )


            # ------------------------------------------------
            # Y axis
            # ------------------------------------------------

            max_count = float(
                window[
                    "count"
                ].max()
            )


            if max_count <= 0:

                max_count = 1


            ax.set_ylim(
                0,
                max_count * 1.1
            )


            fig.tight_layout()


            # ------------------------------------------------
            # Sanitize reference name for filename
            # ------------------------------------------------

            safe_chrom = re.sub(
                r"[^A-Za-z0-9_.-]+",
                "_",
                str(chrom)
            )


            out_pdf = (
                output_dir
                / (
                    f"read_start_"
                    f"{prefix}_"
                    f"{safe_chrom}.pdf"
                )
            )


            fig.savefig(
                out_pdf,
                format="pdf",
                dpi=300,
                bbox_inches="tight"
            )


            plt.close(
                fig
            )


            plots_created += 1


            print(
                f"  Created: "
                f"{out_pdf.name}"
            )


    # ========================================================
    # Finish
    # ========================================================

    print()
    print("========================================")
    print("TSS plotting completed")
    print("----------------------------------------")

    print(
        f"Plots created: "
        f"{plots_created}"
    )

    print()

    print(
        f"Output directory:\n"
        f"{output_dir}"
    )

    print("========================================")


if __name__ == "__main__":
    main()
