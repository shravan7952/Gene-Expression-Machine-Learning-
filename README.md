# Cancer Gene Analysis and Classification

## Overview
This project involves analyzing gene data across various cancer types, visualizing the distribution of specific, overlapping, and common genes, and building a machine learning model to classify cancer types based on RNA data.

---

## Table of Contents
1. [Project Structure](#project-structure)
2. [Data Processing](#data-processing)
3. [Visualization](#visualization)
4. [Machine Learning](#machine-learning)
5. [Feature Importance](#feature-importance)
6. [Dimensionality Reduction (PCA)](#dimensionality-reduction-pca)
7. [Dependencies](#dependencies)

---

## Project Structure

- **`gene_results/`**: Contains text files with specific, overlapping, and common gene data for each cancer type.
- **`plots/`**: Stores output plots generated during analysis.
- **Cancer Data Files**: `.tsv` files with RNA data for various cancers, used for classification.

---

## Data Processing

### Gene Analysis
1. **Common Genes**: Read from `common_genes.txt`.
2. **Specific Genes**: Extracted from files ending with `_specific_genes.txt`.
3. **Overlapping Genes**: Extracted from files ending with `_overlapping_genes.txt`.

### Data Preparation
- Gene counts (specific, overlapping, and common) are computed and aligned across cancer types for visualization.

---

## Visualization

### Stacked Bar Chart
Displays the distribution of common, specific, and overlapping genes for each cancer type.

### Venn Diagram
Shows the overlap between common, specific, and overlapping genes across all cancer types.

### Bar Plots
1. **Specific Genes Per Cancer**: A bar plot for the count of specific genes for each cancer.
2. **Overlapping Genes Analysis**: Visualizes the number of overlapping genes for each cancer type.

### Heatmap
Displays a pairwise overlap matrix of genes between cancer types.

---

## Machine Learning

### Data Preparation
- RNA data from `.tsv` files is merged into a single dataset.
- Cancer type labels are encoded using `LabelEncoder`.
- Data is split into training and testing sets.

### Classification
- A **Random Forest Classifier** is trained on the data to predict cancer types.
- Model performance is evaluated using:
  - **Accuracy**
  - **Classification Report**
  - **Confusion Matrix**

### Model Saving
- The trained model and label encoder are saved using `joblib`.

---

## Feature Importance
A bar chart displays the top 10 most important RNA genes contributing to the cancer type classification.

---

## Dimensionality Reduction (PCA)

### Standardization
Data is scaled using `StandardScaler` for better PCA performance.

### PCA
- Reduces high-dimensional RNA data to 2 principal components.
- A scatter plot visualizes cancer types in the reduced feature space.

---

## Dependencies

The following Python libraries are used in this project:
- **`os`**: File operations.
- **`pandas`**: Data manipulation and analysis.
- **`matplotlib`**: Plotting graphs.
- **`seaborn`**: Statistical data visualization.
- **`matplotlib_venn`**: Creating Venn diagrams.
- **`sklearn`**: Machine learning tools for preprocessing, classification, and evaluation.
- **`joblib`**: Saving and loading models.

---

## Outputs

- **Plots**:
  - Gene distribution bar charts.
  - Venn diagrams for gene category overlap.
  - Heatmap for pairwise gene overlap.
  - Feature importance chart.
  - PCA scatter plot.
- **Model**: Saved Random Forest Classifier and label encoder.

---

## How to Run

1. Place the gene files in the `gene_results` directory.
2. Add `.tsv` RNA data files for cancers in the project directory.
3. Run the script to process data, generate visualizations, and train the model.
4. Check the `plots/` directory for visual outputs.

---

## Notes
- Ensure data files are correctly formatted and complete.
- Update paths and configurations as needed to match your system.

---
