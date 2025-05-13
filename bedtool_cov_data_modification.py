import pandas as pd

# Load the CSV file
file_path = './chromosome_coverage_bedtools.csv'
df = pd.read_csv(file_path)

# Step 1: Remove rows where "Chromosome" contains "genome"
df_filtered = df[~df['Chromosome'].str.contains("genome", case=False, na=False)].copy()

# Step 2: Calculate the median "Mean Coverage" for each "File Name" and normalize "Mean Coverage"
df_filtered['Median Coverage'] = df_filtered.groupby('File Name')['Mean Coverage'].transform('median')
df_filtered['Ploidy'] = (df_filtered['Mean Coverage'] / df_filtered['Median Coverage']) * 2

# Step 3: Round "Ploidy" to two decimal places
df_filtered['Ploidy'] = df_filtered['Ploidy'].round(2)

# Step 4: Rename "File Name" column to "Sample"
df_filtered.rename(columns={'File Name': 'Sample'}, inplace=True)

# Step 5: Transform "Sample" column to keep only the part before "."
df_filtered['Sample'] = df_filtered['Sample'].str.split('.', expand=True)[0]

# Drop the "Median Coverage" column as it was intermediate
df_final = df_filtered.drop(columns=['Median Coverage'])

# Save the final DataFrame to a new CSV file
output_file_path = './processed_chromosome_coverage_with_sample.csv'
df_final.to_csv(output_file_path, index=False)

print("Processed file saved to:", output_file_path)