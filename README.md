# SeuratVizHelper

**SeuratVizHelper** makes single-cell visualization in Seurat easier and more publication-ready. It bridges the gap between R/Seurat and Python/Scanpy aesthetics.

## Features

- 📊 **Stacked Violin Plots**: Scanpy-style publication-ready visualizations
- 🌳 **Automatic Clustering**: Dendrogram-based cluster ordering
- 🎨 **Customizable**: Colors, sizes, and cluster orders
- 💾 **Easy Export**: Built-in high-resolution image saving

## AI Usage

This package was made by AI coding tools (Gemini and Github copilot)

## Installation

Install directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("AyumuOkumura/SeuratVizHelper")
```

## Quick Start

```r
library(Seurat)
library(SeuratVizHelper)

# Load your Seurat object
seurat_obj <- readRDS("seurat_obj.rds")

# Create a stacked violin plot (uses DefaultAssay automatically)
# Save to file
ndim = 20
StackVln(
  seurat_object = seurat_obj,
  features = c("IL7R", "CCR7", "S100A4"),
  group.by = "seurat_clusters",
  dendrogram_method = "dims",
  ndim = ndim,
  plot_width = 12,
  plot_heights = c(0.3, 3),
  # Define the height ratio between the "top graph (dendrogram)" and "bottom graph (violin plot)"
  save_dir = "./figures",
         )

genes <- c("CD3D", "CD8A", "CD4", "MS4A1", "CD14")

StackVln(
  seurat_object = seurat_obj,
  features = genes,
  group.by = "seurat_clusters",
  dendrogram_method = "dims",
  ndim = 30,
  plot_width = 5,
  plot_heights = c(0.3, 3), 
  color_high = "#BD2130"
  save_dir = file.path("png/ViolinPlot")
         )
```

## Advanced Usage
### Custom Colors

```r
StackVln(seurat_obj,
         features = c("CD3D", "CD8A"),
         color_low = "lightblue",
         color_high = "darkred")
```

### Dendrogram Calculation Methods

The `StackVln()` function uses Seurat's `BuildClusterTree()` to calculate hierarchical relationships between clusters. Three methods are available:

#### Method 1: Using Dimensionality Reduction (Recommended)

This method uses PCA or other dimensionality reduction results to calculate cluster relationships. It's computationally efficient and robust to noise.

```r
# Use first 50 PCA dimensions
StackVln(seurat_obj,
         features = c("CD3D", "CD8A", "CD4"),
         dendrogram_method = "dims",
         ndim = 50,
         reduction_for_tree = "pca")

# Use default 30 dimensions (ndim not specified)
StackVln(seurat_obj,
         features = c("CD3D", "CD8A", "CD4"),
         dendrogram_method = "dims")
```

**When to use:**
- Standard workflow for most analyses
- When you have performed PCA or other dimensionality reduction
- For robust clustering relationships across the full transcriptome

**Note:** If `ndim` is not specified, the function defaults to using dimensions 1:30 and displays a message.

#### Method 2: Using Plot Features Only

This method calculates the dendrogram based only on the genes you're plotting.

```r
StackVln(seurat_obj,
         features = c("CD3D", "CD8A", "CD4", "MS4A1", "CD14"),
         dendrogram_method = "features")
```

**When to use:**
- When you want cluster relationships specific to your marker genes
- For focused analysis on a particular gene set
- When your features are carefully selected markers

#### Method 3: Using All Variable Features

This method uses all variable features identified in your Seurat object.

```r
StackVln(seurat_obj,
         features = c("CD3D", "CD8A"),
         dendrogram_method = "all_variable")
```

**When to use:**
- When you want comprehensive gene-based clustering
- As an alternative to dimensionality reduction method
- When you don't have dimensionality reduction results

### Plot Size Customization

Control the dimensions of your plots for optimal visualization and publication quality.

#### Automatic Height Calculation

By default, plot height is automatically calculated based on the number of features:
- Formula: `calc_height = max(5, length(features) * 0.4 + 2)`
- Minimum: 5 inches
- Scales: 0.4 inches per feature + 2 inches base

Examples:
- 3 features → 5 inches (minimum)
- 10 features → 6 inches
- 20 features → 10 inches

#### Adjusting Plot Width

```r
# Wider plot for many clusters
StackVln(seurat_obj,
         features = c("CD3D", "CD8A"),
         plot_width = 15)
```

#### Controlling Dendrogram vs Violin Ratio

The `plot_heights` parameter controls the relative space allocated to dendrogram and violin plots:

```r
# Default: minimal dendrogram, focus on violins
StackVln(seurat_obj, features = genes, plot_heights = c(1, 9))

# Emphasize dendrogram
StackVln(seurat_obj, features = genes, plot_heights = c(2, 8))

# Larger dendrogram for detailed branch structure
StackVln(seurat_obj, features = genes, plot_heights = c(3, 7))

# Minimal dendrogram
StackVln(seurat_obj, features = genes, plot_heights = c(0.5, 9.5))
```

## Function Arguments

- `seurat_object`: Seurat object containing the assay data
- `features`: Character vector of gene names to plot
- `group.by`: Column name in meta.data for grouping (default: "seurat_clusters")
- `cluster_order`: Optional character vector to specify cluster order manually
- `assay`: Assay to use (default: NULL, automatically uses DefaultAssay())
- `color_low`: Color for low expression (default: "white")
- `color_high`: Color for high expression (default: "#BD2130")
- `plot_heights`: Numeric vector (length=2) for relative heights of dendrogram and violin (default: c(1, 9))
- `plot_width`: Width of the saved plot in inches (default: 10)
- `save_path`: Full path to save the file. Overrides save_dir
- `save_dir`: Directory to save the file using a default filename
- `ndim`: Number of dimensions for dendrogram (when method = "dims", default: NULL = 30)
- `dendrogram_method`: Calculation method: "features", "dims" (recommended), or "all_variable" (default: "features")
- `reduction_for_tree`: Reduction to use for "dims" method (default: "pca")

## Citation

If you use SeuratVizHelper in your research, please cite:

```
Okumura, A. (2026). SeuratVizHelper: Enhanced Visualization for Seurat.
R package version 0.0.0.9000.
https://github.com/AyumuOkumura/SeuratVizHelper
```

## License

MIT License

## Issues & Contributions

Please report issues at: https://github.com/AyumuOkumura/SeuratVizHelper/issues
Contributions are welcome via pull requests!

## acknowledgement

- seurat <https://satijalab.org/seurat/>
- scannpy stacked_violin <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.stacked_violin.html>
