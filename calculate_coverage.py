import os
import pysam
import pandas as pd
import argparse

def calculate_mean_coverage(bam_path):
    """
    Calculate mean coverage for each chromosome in a BAM file.

    Parameters:
        bam_path (str): Path to the BAM file.

    Returns:
        dict: A dictionary with chromosomes as keys and mean coverage as values.
    """
    coverage_dict = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for ref in bam.references:
            ref_length = bam.get_reference_length(ref)
            if ref_length == 0:
                coverage_dict[ref] = 0
                continue

            # Get coverage using count_coverage which returns a tuple for A, C, G, T
            coverage = bam.count_coverage(ref)
            # Sum coverage across all bases (A, C, G, T)
            total_coverage_per_base = [sum(bases) for bases in zip(*coverage)]
            # Calculate mean coverage
            mean_cov = sum(total_coverage_per_base) / ref_length
            coverage_dict[ref] = mean_cov

    return coverage_dict

def process_bam_files(directory):
    """
    Process all BAM files in a directory to calculate coverage.

    Parameters:
        directory (str): Path to the directory containing BAM files.

    Returns:
        list: A list of dictionaries containing BAM file name, chromosome, and coverage.
    """
    results = []
    bam_files = [f for f in os.listdir(directory) if f.endswith(".bam")]

    if not bam_files:
        print("No BAM files found in the specified directory.")
        return results

    for bam_file in bam_files:
        bam_path = os.path.join(directory, bam_file)
        print(f"Processing {bam_file}...")
        coverage_dict = calculate_mean_coverage(bam_path)
        for chrom, cov in coverage_dict.items():
            results.append({
                'BAM_File': bam_file,
                'Chromosome': chrom,
                'Mean_Coverage': cov
            })

    return results

def save_to_csv(data, output_path):
    """
    Save coverage data to a CSV file.

    Parameters:
        data (list): List of dictionaries containing coverage data.
        output_path (str): Path to the output CSV file.
    """
    if not data:
        print("No data to save.")
        return

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Coverage data saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Calculate chromosome-wise coverage for BAM files and save to CSV.")
    parser.add_argument('-d', '--directory', required=True, help='Path to the directory containing BAM files.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file.')

    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"The directory {args.directory} does not exist.")
        return

    coverage_data = process_bam_files(args.directory)
    save_to_csv(coverage_data, args.output)

if __name__ == "__main__":
    main()

