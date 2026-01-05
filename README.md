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
  ndim = ndim, # default 30 dimensions (ndim not specified)
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
## Output
{seurat_object}_stack_vln.png

<img src="example/pbmc_stack_vln.png" width="600" alt="PBMC Stacked Violin Plot">

## How StackVln Works: Step-by-Step Computation Details

Understanding the computational workflow behind `StackVln()` helps you make informed decisions about parameters and interpret your results. Here's what happens under the hood:

### Step 1: Data Preparation and Validation

1. **Assay Selection**: If no assay is specified, the function uses `DefaultAssay()` to select the assay
2. **Feature Validation**: The function checks that requested features exist in the selected assay
3. **Group Validation**: Verifies that the `group.by` column exists in the Seurat object's metadata

### Step 2: Cluster Tree Construction

The dendrogram is built using Seurat's `BuildClusterTree()` function. The computation method depends on your `dendrogram_method` parameter:

#### Method A: "dims" (Recommended)
```
1. Extract dimensionality reduction data (e.g., PCA)
   → Uses the specified reduction (default: "pca")
   
2. Calculate cluster centroids in reduced space
   → For each cluster, computes mean coordinates across specified dimensions (default: 1:30)
   
3. Compute pairwise distances between cluster centroids
   → Uses Euclidean distance in the reduced dimensional space
   
4. Perform hierarchical clustering
   → Method: Ward's linkage
   → Creates dendrogram based on cluster similarities
```

**Advantages**: Computationally efficient, robust to noise, captures global transcriptomic relationships

#### Method B: "features"
```
1. Extract expression data for specified features only
   → Uses only the genes you're plotting
   
2. Calculate average expression per cluster
   → For each cluster, computes mean expression of each feature
   
3. Compute correlation/distance between clusters
   → Based solely on the selected marker genes
   
4. Build hierarchical tree
   → Ward's linkage clustering on feature-specific relationships
```

**Advantages**: Shows relationships specific to your markers, useful for focused analysis

#### Method C: "all_variable"
```
1. Identify all variable features in the assay
   → Uses features marked as variable in Seurat object
   
2. Extract expression matrix for all variable features
   → Typically 2000-3000 genes
   
3. Calculate cluster-wise average expression
   → Mean expression across all variable features
   
4. Perform hierarchical clustering
   → Creates comprehensive gene-based dendrogram
```

**Advantages**: Comprehensive gene-based approach, no dimensionality reduction needed

### Step 3: Cluster Ordering

1. **Extract dendrogram structure** from the built cluster tree
2. **Determine optimal leaf order** using hierarchical clustering results
3. **Override with manual order** if `cluster_order` is specified by user

### Step 4: Expression Data Extraction and Scaling

```
For each feature:
1. Extract raw expression values from the specified assay
2. Group cells by cluster (using group.by parameter)
3. Scale expression to 0-1 range within each feature
   → Formula: (value - min) / (max - min)
   → This enables color mapping across different expression scales
```

### Step 5: Violin Plot Generation

```
For each feature × cluster combination:
1. Extract expression distribution for all cells in that cluster
2. Calculate kernel density estimation
   → Creates smooth distribution curve
3. Mirror the density plot to create violin shape
4. Color violin by mean expression (scaled 0-1)
   → Low expression → color_low (default: white)
   → High expression → color_high (default: #BD2130)
```

### Step 6: Layout Assembly

1. **Create dendrogram plot** with height ratio from `plot_heights[1]`
2. **Create stacked violin plots** with height ratio from `plot_heights[2]`
3. **Align plots** vertically using `patchwork` package
4. **Calculate total plot height**: 
   - If not specified: `max(5, length(features) × 0.4 + 2)` inches
   - Ensures adequate space for all features

### Step 7: Export

1. **Create output directory** if it doesn't exist
2. **Save plot** as PNG with specified dimensions
   - Default filename: `{seurat_object_name}_stack_vln.png`
   - Resolution: High-quality for publication

### Key Computational Notes

- **Memory efficiency**: The function processes one feature at a time rather than loading all data at once
- **Scaling strategy**: Per-feature scaling (not global) ensures each gene is visually comparable
- **Clustering algorithm**: Uses Ward's criterion (`ward.D2`) which minimizes within-cluster variance
- **Distance metric**: Euclidean distance in PCA space (for "dims" method) or expression space (other methods)

## Advanced Usage
### Custom Colors

```r
StackVln(
  seurat_object = seurat_obj,
  features = c("IL7R", "CCR7", "S100A4"),
  group.by = "seurat_clusters",
  dendrogram_method = "dims",
  plot_width = 12,
  plot_heights = c(0.3, 3),
  color_low = "lightblue",
　color_high = "darkred",
  save_dir = "./figures",
         )
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
dendrogram_method = "features"

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
