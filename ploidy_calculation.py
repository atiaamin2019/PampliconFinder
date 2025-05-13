import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the dataset
csv_file_path = './processed_chromosome_coverage_with_sample.csv'  # Replace with your actual file path
data = pd.read_csv(csv_file_path)

# Pivot the data to create a matrix suitable for a heatmap
heatmap_data = data.pivot(index='Chromosome', columns='Sample', values='Ploidy')

# Rename chromosomes sequentially from Chr1 to Chr35
new_chromosome_names = [f"Chr{i}" for i in range(1, len(heatmap_data) + 1)]
heatmap_data.index = new_chromosome_names

# Sort Sample names to group '_P0' samples together
samples = heatmap_data.columns
p0_samples = [sample for sample in samples if '_P0' in sample]
non_p0_samples = [sample for sample in samples if '_P0' not in sample]

# Rearrange columns with '_P0' samples first, adding an empty column after the first 10 '_P0' samples
split_index = 10 if len(p0_samples) > 10 else len(p0_samples)
new_column_order = p0_samples[:split_index] + ['Empty'] + p0_samples[split_index:] + non_p0_samples
heatmap_data['Empty'] = np.nan  # Add an empty column
heatmap_data = heatmap_data[new_column_order]

# Adjust annotations: show floats and conditionally annotate
annot_data = heatmap_data.applymap(
    lambda x: f"{x:.1f}" if pd.notnull(x) and (x <= 1.5 or x >= 2.5) else ""
)

# Plot the heatmap
plt.figure(figsize=(14, max(10, len(heatmap_data) * 0.5)))

sns.heatmap(
    heatmap_data, 
    annot=annot_data, 
    fmt="",
    cmap="viridis",  
    #cmap="coolwarm", 
    cbar=True, 
    linewidths=0.5,
    linecolor="black",  # Set box border color to black
    mask=heatmap_data.isnull(),  # Mask the empty column for no color/line
)
plt.title('Heatmap of Ploidy Values with Conditional Annotations and Adjusted Columns')
plt.xlabel('Samples')
plt.ylabel('Chromosomes')
plt.xticks(rotation=90)  # Rotate column names to 90 degrees
plt.yticks(rotation=0)
plt.tight_layout()
# Save the figure
output_file_path = './heatmap.png'  
plt.savefig(output_file_path, dpi=300, bbox_inches='tight')

plt.show()


