import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib_venn import venn3

gene_folder = "gene_results"
output_folder = "plots"
os.makedirs(output_folder, exist_ok=True)

# Load common genes
with open(os.path.join(gene_folder, "common_genes.txt")) as f:
    common_genes = set(f.read().splitlines())

# Load specific and overlapping genes
specific_genes = {}
overlapping_genes = {}
for file in os.listdir(gene_folder):
    if file.endswith("_specific_genes.txt"):
        cancer_type = file.split("_specific")[0]
        with open(os.path.join(gene_folder, file)) as f:
            genes = set(f.read().splitlines())
        if "No specific genes for this cancer." not in genes:
            specific_genes[cancer_type] = genes
    elif file.endswith("_overlapping_genes.txt"):
        cancer_type = file.split("_overlapping")[0]
        with open(os.path.join(gene_folder, file)) as f:
            genes = set(f.read().splitlines())
        if "No specific genes for this cancer." in genes:
            genes.remove("No specific genes for this cancer.")
        if genes:
            overlapping_genes[cancer_type] = genes

print(f"Loaded {len(common_genes)} common genes.")
print(f"Loaded specific genes for {len(specific_genes)} cancers.")
print(f"Loaded overlapping genes for {len(overlapping_genes)} cancers.")

# Count specific and overlapping genes
specific_counts = {cancer: len(genes) for cancer, genes in specific_genes.items()}
overlapping_counts = {cancer: len(genes) for cancer, genes in overlapping_genes.items()}

# Merge counts
all_cancer_types = set(specific_counts.keys()).union(set(overlapping_counts.keys()))
aligned_specific_counts = {cancer: specific_counts.get(cancer, 0) for cancer in all_cancer_types}
aligned_overlapping_counts = {cancer: overlapping_counts.get(cancer, 0) for cancer in all_cancer_types}

common_genes_count = len(common_genes)
aligned_common_counts = {cancer: common_genes_count for cancer in all_cancer_types}

sorted_cancers = sorted(all_cancer_types)

# Prepare data for plotting
specific_values = [aligned_specific_counts[cancer] for cancer in sorted_cancers]
overlapping_values = [aligned_overlapping_counts[cancer] for cancer in sorted_cancers]
common_values = [aligned_common_counts[cancer] for cancer in sorted_cancers]

# Plot stacked bar chart
plt.figure(figsize=(12, 12))
plt.barh(sorted_cancers, specific_values, label="Specific Genes", alpha=0.7)
plt.barh(sorted_cancers, overlapping_values, label="Overlapping Genes", alpha=0.7, left=specific_values)
plt.barh(sorted_cancers, common_values, label="Common Genes", alpha=0.7, left=[s + o for s, o in zip(specific_values, overlapping_values)])

plt.xlabel("Number of Genes")
plt.ylabel("Cancer Types")
plt.title("Distribution of Genes Across Cancers")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "gene_distribution_plot.png"), dpi=300)
plt.show()

# Venn Diagram
venn3([set(common_genes), set.union(*specific_genes.values()), set.union(*overlapping_genes.values())],
      ('Common Genes', 'Specific Genes', 'Overlapping Genes'))
plt.title("Gene Category Overlap")
plt.savefig(os.path.join(output_folder, "common_specific_overlapping_genes.png"), dpi=300)
plt.show()

# Bar plot for specific genes
specific_df = pd.DataFrame.from_dict(specific_counts, orient='index', columns=['Specific Gene Count'])
specific_df = specific_df.sort_values(by='Specific Gene Count', ascending=False)

plt.figure(figsize=(12, 6))
sns.barplot(x=specific_df.index, y=specific_df['Specific Gene Count'], edgecolor="black", color="lightblue")
plt.title("Number of Specific Genes for Each Cancer Type")
plt.xlabel("Cancer Types")
plt.ylabel("Number of Specific Genes")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "specific_genes_per_cancer.png"), dpi=300)
plt.show()

# Overlapping genes analysis
overlap_counts = {}
for cancer_type, genes in overlapping_genes.items():
    other_cancers = set(overlapping_genes.keys()) - {cancer_type}
    overlap_counts[cancer_type] = sum(len(genes.intersection(overlapping_genes[other])) for other in other_cancers)

overlap_df = pd.DataFrame.from_dict(overlap_counts, orient='index', columns=['Overlap Count'])
overlap_df = overlap_df.sort_values(by='Overlap Count', ascending=False)

# Bar plot for overlapping genes
# Create a bar plot for overlapping genes
plt.figure(figsize=(18, 10))
sns.barplot(
    x=overlap_df.index,
    y=overlap_df['Overlap Count'],
    palette=sns.color_palette("viridis", len(overlap_df)),
    edgecolor="black"
)
plt.title("Number of Overlapping Genes for Each Cancer Type")
plt.xlabel("Cancer Types")
plt.ylabel("Number of Overlapping Genes")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Construct the output path
output_path = os.path.join(output_folder, "overlapping_genes_all_cancers_barplot.png")
print(f"Saving plot to: {output_path}")  # Debugging line to check the output path

# Save the plot
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

# Pairwise overlap matrix
cancer_types = list(overlapping_genes.keys())
overlap_matrix = pd.DataFrame(index=cancer_types, columns=cancer_types, dtype=int)

for cancer1 in cancer_types:
    for cancer2 in cancer_types:
        overlap_matrix.loc[cancer1, cancer2] = len(overlapping_genes[cancer1].intersection(overlapping_genes[cancer2]))

overlap_matrix.fillna(0, inplace=True)
overlap_matrix = overlap_matrix.astype(int)

# Heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(overlap_matrix, annot=True, cmap="coolwarm", fmt="d", cbar_kws={'label': 'Number of Overlapping Genes'})
plt.title("Pairwise Overlap of Genes Between Cancer Types")
plt.xlabel("Cancer Types")
plt.ylabel("Cancer Types")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "pairwise_overlap_heatmap_fixed.png"), dpi=300)
plt.show()