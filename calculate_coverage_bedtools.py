import os
import subprocess
import pandas as pd

# Directory containing BAM files
bam_directory = "./"
output_csv = "chromosome_coverage_bedtools.csv"
genome_file = "/Users/aamin/Desktop/Leishmania_Project/Amplicon_Project/L.braziliensis/Reference_Genome/ncbi_dataset_M2904/ncbi_dataset/data/GCA_900537975.2/GCA_900537975.2_LBRM_annotationDEFINITIVO_genomic.fna.fai"  # Reference genome file (.fai)

# List to store the results
coverage_data = []

# Iterate over all BAM files in the directory
for bam_file in os.listdir(bam_directory):
    if bam_file.endswith(".bam"):
        bam_path = os.path.join(bam_directory, bam_file)
        
        # Command to calculate coverage using bedtools genomecov
        bedtools_command = f"bedtools genomecov -ibam {bam_path} -g {genome_file}"
        result = subprocess.run(bedtools_command, shell=True, capture_output=True, text=True)

        # Check if the command was successful
        if result.returncode != 0:
            print(f"Error processing {bam_file}: {result.stderr}")
            continue
        
        # Parse the bedtools output
        coverage_by_chromosome = {}
        for line in result.stdout.splitlines():
            fields = line.split()
            chromosome = fields[0]
            depth = int(fields[1])
            num_bases = int(fields[2])
            if chromosome not in coverage_by_chromosome:
                coverage_by_chromosome[chromosome] = {"total_coverage": 0, "total_bases": 0}
            coverage_by_chromosome[chromosome]["total_coverage"] += depth * num_bases
            coverage_by_chromosome[chromosome]["total_bases"] += num_bases

        # Calculate mean coverage for each chromosome
        for chromosome, stats in coverage_by_chromosome.items():
            if stats["total_bases"] > 0:
                mean_coverage = stats["total_coverage"] / stats["total_bases"]
                coverage_data.append([bam_file, chromosome, mean_coverage])

# Create a DataFrame and save to CSV
coverage_df = pd.DataFrame(coverage_data, columns=["File Name", "Chromosome", "Mean Coverage"])
coverage_df.to_csv(output_csv, index=False)

print(f"Coverage data saved to {output_csv}")