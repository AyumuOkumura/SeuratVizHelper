# StackVln 0.0.1.9000

## Major Changes

* **BREAKING**: Replaced `hclust()` with Seurat's `BuildClusterTree()` for dendrogram calculation
* Added `dendrogram_method` parameter with three options: "features", "dims", "all_variable"
* Added `ndim` parameter to specify number of dimensions for PCA-based dendrogram
* Added `reduction_for_tree` parameter to specify reduction method (default: "pca")

## Improvements

* More consistent with Seurat's standard workflow
* Better integration with Seurat's dimensional reduction results
* Informative message when using default ndim = 30
* Enhanced README with comprehensive documentation on:
  - Dendrogram calculation methods and when to use each
  - Plot size customization for publication-quality figures
  - Best practices for journal submission

## Documentation

* Added detailed explanation of dendrogram calculation methods
* Added plot size customization guide
* Added publication-quality figure examples

---

# StackVln 0.0.0.9000

## New Features

* Initial release
* `StackVln()` function for Scanpy-style stacked violin plots
* Automatic dendrogram-based cluster ordering
* Flexible color customization
* Built-in high-resolution plot saving

## Improvements

* Comprehensive input validation
* Automatic directory creation for saved plots
* Enhanced error handling and informative messages
* Detailed documentation with usage examples

## Dependencies

* Seurat >= 4.0.0
* ggplot2, dplyr, tidyr, patchwork, ggdendro, tibble
* rlang, magrittr
