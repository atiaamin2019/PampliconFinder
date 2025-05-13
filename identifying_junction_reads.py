import re
import os
import gzip

# Define the query sequences
#query sequences qPCR primers for LINF  mrpa
query1_sequence = "AGGAAAACAATGAAAAGATAG"
query2_sequence = "CGCGGCCGACAGAATAATAA"

query3_sequence = "CTATCTTTTCATTGTTTTCCT"
query4_sequence = "TTATTATTCTGTCGGCCGCG"

#for linear amplicon searching switching the directions of the probes
#query1_sequence = "AGGAAAACAATGAAAAGATAG"
#query4_sequence = "CGCGGCCGACAGAATAATAA"

#query3_sequence = "CTATCTTTTCATTGTTTTCCT"
#query2_sequence = "TTATTATTCTGTCGGCCGCG"



# Create patterns for both sets of sequences
pattern1 = re.compile(
    re.escape(query1_sequence) + ".*" + re.escape(query2_sequence) +
    "|" +
    re.escape(query2_sequence) + ".*" + re.escape(query1_sequence)
)
pattern2 = re.compile(
    re.escape(query3_sequence) + ".*" + re.escape(query4_sequence) +
    "|" +
    re.escape(query4_sequence) + ".*" + re.escape(query3_sequence)
)

# Directory containing the FASTQ.GZ files
input_dir = "./"  # Change to your actual directory

# Process each FASTQ.GZ file in the directory
for file in os.listdir(input_dir):
    if file.endswith(".fastq.gz"):
        input_file = os.path.join(input_dir, file)
        output_file = os.path.join(input_dir, f"{file.rsplit('.fastq.gz', 1)[0]}.mrpA_qpcr_product_same_direction.fasta")

        with gzip.open(input_file, "rt") as input_handle, open(output_file, "w") as output_handle:
            while True:
                # Read each set of 4 lines (one FASTQ record)
                header = input_handle.readline().strip()
                sequence = input_handle.readline().strip()
                plus = input_handle.readline().strip()
                quality = input_handle.readline().strip()

                # Check if we reached the end of the file
                if not header:
                    break

                # Search for patterns in the sequence
                match1 = pattern1.search(sequence)
                match2 = pattern2.search(sequence)

                # Write the matching sequence in FASTA format
                if match1:
                    trimmed_sequence = match1.group(0)
                    output_handle.write(f">{header[1:]}\n{trimmed_sequence}\n")
                elif match2:
                    trimmed_sequence = match2.group(0)
                    output_handle.write(f">{header[1:]}\n{trimmed_sequence}\n")

        print(f"Processed {file} and saved results to {output_file}")

