import gzip
import os
from glob import glob

def parse_txt_file(txt_file_path):
    """
    Parses the TXT file to extract names where length > 300.
    """
    filtered_names = set()
    with open(txt_file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip the header if present
                continue
            name, length = line.strip().split('\t')
            if int(length) > 300:
                filtered_names.add(name)
    return filtered_names

def extract_sequences_from_fastq(fastq_file_path, filtered_names, output_fasta_path):
    """
    Extracts sequences from a fastq.gz file that match filtered_names and saves them to a fasta file.
    """
    with gzip.open(fastq_file_path, 'rt') as fastq, open(output_fasta_path, 'w') as fasta:
        write_seq = False
        for line in fastq:
            line = line.strip()
            if line.startswith('@'):  # Header line in fastq
                seq_name = line[1:]  # Remove '@' from the header
                if seq_name in filtered_names:
                    write_seq = True
                    fasta.write(f">{seq_name}\n")  # Write header in fasta format
                else:
                    write_seq = False
            elif write_seq:  # Sequence line
                fasta.write(f"{line}\n")


def process_files(input_dir, output_dir):
    """
    Processes all .txt and .fastq.gz files in the input directory and saves results in the output directory.
    """
    # Find all .txt files in the input directory
    txt_files = glob(os.path.join(input_dir, "*.txt"))
    
    for txt_file in txt_files:
        # Infer the corresponding fastq.gz file from the txt file name
        base_name = os.path.splitext(os.path.basename(txt_file))[0]
        fastq_file = os.path.join(input_dir, f"{base_name}.fastq.gz")
        
        if not os.path.exists(fastq_file):
            print(f"FASTQ file not found for {txt_file}. Skipping...")
            continue
        
        # Output FASTA file path
        output_fasta_file = os.path.join(output_dir, f"{base_name}_filtered_sequences.fasta")
        
        print(f"Processing: {txt_file} and {fastq_file}")
        
        # Step 1: Parse the TXT file to get filtered names
        filtered_names = parse_txt_file(txt_file)
        
        # Step 2: Extract sequences from the fastq.gz file
        extract_sequences_from_fastq(fastq_file, filtered_names, output_fasta_file)
        
        print(f"Filtered sequences saved to {output_fasta_file}")

# Example usage:
input_directory = "./"  # Replace with your input directory path
output_directory = "./extracted_reads/"  # Replace with your output directory path

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Process all files in the input directory
process_files(input_directory, output_directory)
