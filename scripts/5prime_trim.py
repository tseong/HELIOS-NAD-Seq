import regex as re
import time
import os
import multiprocessing
from functools import partial
import glob

# Start the timer
start_time = time.perf_counter()

def trim(input_file, output_files, adapter_keys=None):
    total_reads = 0
    short_reads = 0
    forward_trimmed_5prime = 0
    forward_tot_length = 0

    bc, file_names = input_file

    # Dictionary of adapters (only forward_5prime)
    adapters = {
        'forward_5prime': '^[GN]*',
    }

    # Initialize compiled adapters with default empty strings
    compiled_adapters = {key: re.compile('') for key in adapters}

    if adapter_keys:
        if isinstance(adapter_keys, str):
            adapter_keys = [adapter_keys]  # Convert to list for uniform processing

        for key in adapter_keys:
            if key in adapters:
                compiled_adapters[key] = re.compile(adapters[key])

    with open(file_names['forward'], 'r') as forward_file, open(file_names['reverse'], 'r') as reverse_file:
        buffer_size = 25000  # Number of records to accumulate before writing
        forward_buffer = []
        reverse_buffer = []

        while True:
            forward_lines = [next(forward_file, None) for _ in range(4)]
            reverse_lines = [next(reverse_file, None) for _ in range(4)]

            if None in forward_lines or None in reverse_lines:
                break

            total_reads += 1
            forward_seq = forward_lines[1].strip()
            trimmed_5prime = False
            adapter_object = compiled_adapters.get('forward_5prime')

            if adapter_object.pattern == '':
                forward_lines[1] = forward_lines[1]
                forward_lines[3] = forward_lines[3]
            else:
                match_forward_5prime = compiled_adapters.get('forward_5prime').search(forward_lines[1])
                if match_forward_5prime:
                    end_index = match_forward_5prime.end()
                    forward_lines[1] = forward_lines[1][end_index:]
                    forward_lines[3] = forward_lines[3][end_index:]
                    if end_index != 0:
                        forward_trimmed_5prime += 1
                        trimmed_5prime = True

            if len(forward_lines[1].strip()) < 18:
                short_reads += 1
                          forward_tot_length += len(forward_lines[1].strip())

            forward_buffer.extend(forward_lines)
            reverse_buffer.extend(reverse_lines)

            if len(forward_buffer) >= buffer_size * 4:
                with open(output_files[bc]['forward'], 'a') as forward_out, open(output_files[bc]['reverse'], 'a') as reverse_out:
                    forward_out.writelines(forward_buffer)
                    reverse_out.writelines(reverse_buffer)
                    forward_buffer = []
                    reverse_buffer = []

        if forward_buffer:
            with open(output_files[bc]['forward'], 'a') as forward_out, open(output_files[bc]['reverse'], 'a') as reverse_out:
                forward_out.writelines(forward_buffer)
                reverse_out.writelines(reverse_buffer)

    # Calculate and print percentages
    forward_average_length = forward_tot_length / total_reads
    percent_forward_trimmed_5prime = (forward_trimmed_5prime / total_reads) * 100 if total_reads > 0 else 0
    percent_short = (short_reads / total_reads) * 100 if total_reads > 0 else 0

    print(f"{bc} forward insert average length: {forward_average_length}")
    print(f"{bc} percentage of forward reads trimmed by 5 prime adapter: {percent_forward_trimmed_5prime}%")
    print(f"{bc} percentage short (shorter than 18 bp): {percent_short}%")

# Path to your specified folder
folder_path = '/gpfs/bwfor/work/ws/hd_uv268-YZ817_eColiHelios_2/cutadapt/'
all_forward_files = glob.glob(f'{folder_path}*R1*')
all_reverse_files = glob.glob(f'{folder_path}*R2*')

# Generate input_files dictionary based on the files found
input_files = {
    f.split('/')[-1]: {
        'forward': f,
        'reverse': f.replace('R1_001.fastq', 'R2_001.fastq')
    }
    for f in all_forward_files if f.replace('R1_001.fastq', 'R2_001.fastq') in all_reverse_files
}

# Assuming you want to append '.custom_trimmed' before the file extension for output files
output_files = {
    name: {
        'forward': f['forward'].replace('R1_001.fastq', 'R1_001_trimmed.fastq'),
        'reverse': f['reverse'].replace('R2_001.fastq', 'R2_001_trimmed.fastq')
    }
    for name, f in input_files.items()
}

# Create a function with the output_files and adapter_keys argument prefilled because pool.map can only take one argument from trim function
trim_partial = partial(trim, output_files=output_files, adapter_keys=['forward_5prime'])

def main():
    # Place the main script logic here
    # Generate partial function and run multiprocessing pool
    trim_partial = partial(trim, output_files=output_files, adapter_keys=['forward_5prime'])

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(trim_partial, input_files.items())

    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"The code took {elapsed_time} seconds to execute.")

if __name__ == "__main__":
    main()
