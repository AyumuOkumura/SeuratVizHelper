# SeuratVizHelper

**SeuratVizHelper** makes single-cell visualization in Seurat easier and more publication-ready. It bridges the gap between R/Seurat and Python/Scanpy aesthetics.

## Features

- 📊 **Stacked Violin Plots**: Scanpy-style publication-ready visualizations
- 🌳 **Automatic Clustering**: Dendrogram-based cluster ordering
- 🎨 **Customizable**: Colors, sizes, and cluster orders
- 💾 **Easy Export**: Built-in high-resolution image saving

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
pbmc <- readRDS("pbmc.rds")

# Create a stacked violin plot
StackVln(pbmc, 
         features = c("CD3D", "CD8A", "CD4", "MS4A1", "CD14"),
         assay = "RNA",
         color_high = "#BD2130")

# Save to file
StackVln(pbmc, 
         features = c("IL7R", "CCR7", "S100A4"),
         save_dir = "./figures",
         plot_width = 12)
```

## Advanced Usage

### Custom Cluster Order

```r
StackVln(pbmc,
         features = c("CD79A", "MS4A1"),
         cluster_order = c("3", "0", "1", "2", "4"))
```

### Different Assay and Grouping

```r
StackVln(pbmc,
         features = c("nFeature_RNA", "nCount_RNA"),
         group.by = "cell_type",
         assay = "SCT")
```

### Custom Colors

```r
StackVln(pbmc,
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
StackVln(pbmc,
         features = c("CD3D", "CD8A", "CD4"),
         dendrogram_method = "dims",
         ndim = 50,
         reduction_for_tree = "pca")

# Use default 30 dimensions (ndim not specified)
StackVln(pbmc,
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
StackVln(pbmc,
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
StackVln(pbmc,
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
StackVln(pbmc,
         features = c("CD3D", "CD8A"),
         plot_width = 15)
```

#### Controlling Dendrogram vs Violin Ratio

The `plot_heights` parameter controls the relative space allocated to dendrogram and violin plots:

```r
# Default: minimal dendrogram, focus on violins
StackVln(pbmc, features = genes, plot_heights = c(1, 9))

# Emphasize dendrogram
StackVln(pbmc, features = genes, plot_heights = c(2, 8))

# Larger dendrogram for detailed branch structure
StackVln(pbmc, features = genes, plot_heights = c(3, 7))

# Minimal dendrogram
StackVln(pbmc, features = genes, plot_heights = c(0.5, 9.5))
```

#### Publication-Quality Settings

```r
# For journal submission (Nature/Cell style)
StackVln(pbmc,
         features = c("CD3D", "CD8A", "CD4", "MS4A1", "CD14",
                      "FCGR3A", "LYZ", "PPBP", "IL7R", "CCR7"),
         plot_width = 15,
         plot_heights = c(1.5, 8.5),
         dendrogram_method = "dims",
         ndim = 50,
         save_dir = "./figures")
```

**Tips:**
- Check journal requirements for figure dimensions
- Use DPI 300 for publication quality (automatically set in ggsave)
- Wider plots work better for many clusters (>8)
- Taller plots work better for many features (>10)

## Function Arguments

- `seurat_object`: Seurat object containing the assay data
- `features`: Character vector of gene names to plot
- `group.by`: Column name in meta.data for grouping (default: "seurat_clusters")
- `cluster_order`: Optional character vector to specify cluster order manually
- `assay`: Assay to use (default: "SCT")
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

