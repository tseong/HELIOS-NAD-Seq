import regex as re
import time
import os
import multiprocessing
from functools import partial
import glob
import argparse


def trim(input_file, output_files, adapter_keys=None):
    total_reads = 0
    short_reads = 0
    forward_trimmed_5prime = 0
    forward_tot_length = 0

    bc, file_names = input_file

    # Dictionary of adapters
    adapters = {
        "forward_5prime": r"^[GN]*",
    }

    # Initialize compiled adapters
    compiled_adapters = {key: re.compile("") for key in adapters}

    if adapter_keys:
        if isinstance(adapter_keys, str):
            adapter_keys = [adapter_keys]

        for key in adapter_keys:
            if key in adapters:
                compiled_adapters[key] = re.compile(adapters[key])

    forward_buffer = []
    reverse_buffer = []
    buffer_size = 25000

    with open(file_names["forward"], "r") as forward_file, \
         open(file_names["reverse"], "r") as reverse_file:

        while True:
            forward_lines = [next(forward_file, None) for _ in range(4)]
            reverse_lines = [next(reverse_file, None) for _ in range(4)]

            if None in forward_lines or None in reverse_lines:
                break

            total_reads += 1

            adapter_object = compiled_adapters["forward_5prime"]

            if adapter_object.pattern != "":
                match_forward_5prime = adapter_object.search(
                    forward_lines[1].strip()
                )

                if match_forward_5prime:
                    end_index = match_forward_5prime.end()

                    # Trim sequence
                    sequence = forward_lines[1].strip()
                    quality = forward_lines[3].strip()

                    sequence = sequence[end_index:]
                    quality = quality[end_index:]

                    forward_lines[1] = sequence + "\n"
                    forward_lines[3] = quality + "\n"

                    if end_index != 0:
                        forward_trimmed_5prime += 1

            # Length after trimming
            trimmed_length = len(forward_lines[1].strip())
            forward_tot_length += trimmed_length

            if trimmed_length < 18:
                short_reads += 1

            forward_buffer.extend(forward_lines)
            reverse_buffer.extend(reverse_lines)

            if len(forward_buffer) >= buffer_size * 4:
                with open(output_files[bc]["forward"], "a") as forward_out, \
                     open(output_files[bc]["reverse"], "a") as reverse_out:

                    forward_out.writelines(forward_buffer)
                    reverse_out.writelines(reverse_buffer)

                forward_buffer = []
                reverse_buffer = []

        # Write remaining reads
        if forward_buffer:
            with open(output_files[bc]["forward"], "a") as forward_out, \
                 open(output_files[bc]["reverse"], "a") as reverse_out:

                forward_out.writelines(forward_buffer)
                reverse_out.writelines(reverse_buffer)

    # Calculate statistics
    if total_reads > 0:
        forward_average_length = forward_tot_length / total_reads
        percent_forward_trimmed_5prime = (
            forward_trimmed_5prime / total_reads
        ) * 100
        percent_short = (short_reads / total_reads) * 100
    else:
        forward_average_length = 0
        percent_forward_trimmed_5prime = 0
        percent_short = 0

    print(
        f"{bc} forward insert average length: "
        f"{forward_average_length:.2f}"
    )
    print(
        f"{bc} percentage of forward reads trimmed by "
        f"5 prime adapter: {percent_forward_trimmed_5prime:.2f}%"
    )
    print(
        f"{bc} percentage short (shorter than 18 bp): "
        f"{percent_short:.2f}%"
    )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Trim leading G/N bases from the 5' end of "
            "demultiplexed HELIOS FASTQ files."
        )
    )

    parser.add_argument(
        "input_dir",
        help="Directory containing Cutadapt-demultiplexed FASTQ files"
    )

    parser.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help=(
            "Directory for trimmed FASTQ files. "
            "Default: <parent_of_input_dir>/5prime_trimmed"
        )
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="Number of parallel processes"
    )

    args = parser.parse_args()

    start_time = time.perf_counter()

    # Convert to absolute paths
    input_dir = os.path.abspath(args.input_dir)

    if not os.path.isdir(input_dir):
        raise FileNotFoundError(
            f"Input directory does not exist: {input_dir}"
        )

    # Set output directory
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
    else:
        output_dir = os.path.join(
            os.path.dirname(input_dir),
            "5prime_trimmed"
        )

    os.makedirs(output_dir, exist_ok=True)

    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Processes:        {args.threads}")

    # Find R1 and R2 files
    all_forward_files = sorted(
        glob.glob(os.path.join(input_dir, "*R1*.fastq"))
    )

    all_reverse_files = set(
        glob.glob(os.path.join(input_dir, "*R2*.fastq"))
    )

    if not all_forward_files:
        raise FileNotFoundError(
            f"No R1 FASTQ files found in {input_dir}"
        )

    # Pair R1 and R2 files
    input_files = {}

    for forward_file in all_forward_files:

        reverse_file = forward_file.replace(
            "R1.fastq",
            "R2.fastq"
        )

        if reverse_file in all_reverse_files:
            sample_name = os.path.basename(forward_file)

            input_files[sample_name] = {
                "forward": forward_file,
                "reverse": reverse_file
            }
        else:
            print(
                f"WARNING: No matching R2 file found for "
                f"{forward_file}"
            )

    if not input_files:
        raise RuntimeError(
            "No valid R1/R2 FASTQ pairs were found."
        )

    print(f"FASTQ pairs found: {len(input_files)}")

    # Define output files
    output_files = {}

    for name, files in input_files.items():

        forward_basename = os.path.basename(
            files["forward"]
        ).replace(
            "R1.fastq",
            "R1_trimmed.fastq"
        )

        reverse_basename = os.path.basename(
            files["reverse"]
        ).replace(
            "R2.fastq",
            "R2_trimmed.fastq"
        )

        output_files[name] = {
            "forward": os.path.join(
                output_dir,
                forward_basename
            ),
            "reverse": os.path.join(
                output_dir,
                reverse_basename
            )
        }

        # Prevent accidental appending to previous output
        for output_file in output_files[name].values():
            if os.path.exists(output_file):
                os.remove(output_file)

    trim_partial = partial(
        trim,
        output_files=output_files,
        adapter_keys=["forward_5prime"]
    )

    with multiprocessing.Pool(
        processes=args.threads
    ) as pool:

        pool.map(
            trim_partial,
            input_files.items()
        )

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time

    print(
        f"5' trimming completed in "
        f"{elapsed_time:.2f} seconds."
    )


if __name__ == "__main__":
    main()
