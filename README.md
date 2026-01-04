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

