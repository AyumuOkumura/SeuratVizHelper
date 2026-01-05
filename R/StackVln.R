#' Stacked Violin Plot (Scanpy Style)
#'
#' SeuratオブジェクトからScanpy風のスタックバイオリンプロットを描画します。
#' デンドログラムによるクラスタ順序の自動調整と、画像保存機能を備えています。
#'
#' @param seurat_object Seurat object containing the assay data.
#' @param features Vector of gene names to plot.
#' @param group.by Column name in meta.data for grouping (default: "seurat_clusters").
#' @param cluster_order Optional vector to specify cluster order manually.
#' @param assay Assay to use (default: NULL, uses DefaultAssay()).
#' @param color_low Color for low expression (default: "white").
#' @param color_high Color for high expression (default: "#BD2130").
#' @param plot_heights Vector (len=2) relative heights of dendrogram and violin (default: c(1, 9)).
#' @param plot_width Width of the saved plot in inches (default: 10).
#' @param save_path Full path to save the file. Overrides save_dir.
#' @param save_dir Directory to save the file using a default filename.
#' @param ndim Number of dimensions to use when dendrogram_method = "dims". If NULL, uses default 1:30.
#' @param dendrogram_method Method for calculating dendrogram. Options: "features" (uses specified features), 
#'   "dims" (uses dimensionality reduction - recommended), "all_variable" (uses all variable features).
#' @param reduction_for_tree Name of reduction to use when dendrogram_method = "dims" (default: "pca").
#'
#' @return A patchwork object combining the dendrogram and violin plots.
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(SeuratVizHelper)
#' 
#' # Basic usage
#' StackVln(pbmc, features = c("CD3D", "CD8A", "CD4"))
#' 
#' # Custom cluster order and colors
#' StackVln(pbmc, 
#'          features = c("MS4A1", "CD79A"),
#'          cluster_order = c("0", "3", "1", "2"),
#'          color_high = "darkblue")
#' 
#' # Save to file
#' StackVln(pbmc, 
#'          features = c("CD14", "LYZ"),
#'          save_dir = "./plots")
#' }
#'
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import ggdendro
#' @import patchwork
#' @importFrom rlang !! sym
StackVln <- function(
    seurat_object,
    features,
    group.by = "seurat_clusters",
    cluster_order = NULL,
    assay = NULL,
    color_low = "white",
    color_high = "#BD2130", 
    plot_heights = c(1, 9),
    plot_width = 10,
    save_path = NULL,
    save_dir = NULL,
    ndim = NULL,
    dendrogram_method = c("features", "dims", "all_variable"),
    reduction_for_tree = "pca"
) {
  
  # Auto-detect assay if not specified
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_object)
    message(sprintf("Using default assay: %s", assay))
  }
  
  # Input validation
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }
  
  if (!is.character(features) || length(features) == 0) {
    stop("features must be a non-empty character vector")
  }
  
  if (!assay %in% names(seurat_object@assays)) {
    stop(sprintf("Assay '%s' not found in Seurat object. Available assays: %s",
                 assay, paste(names(seurat_object@assays), collapse = ", ")))
  }
  
  if (!group.by %in% colnames(seurat_object@meta.data)) {
    stop(sprintf("Column '%s' not found in meta.data. Available columns: %s",
                 group.by, paste(head(colnames(seurat_object@meta.data), 10), collapse = ", ")))
  }
  
  if (length(plot_heights) != 2 || !is.numeric(plot_heights)) {
    stop("plot_heights must be a numeric vector of length 2")
  }
  
  if (!is.null(cluster_order)) {
    if (!is.character(cluster_order) && !is.numeric(cluster_order)) {
      stop("cluster_order must be a character or numeric vector")
    }
    cluster_order <- as.character(cluster_order)
    available_clusters <- unique(as.character(seurat_object@meta.data[[group.by]]))
    invalid_clusters <- setdiff(cluster_order, available_clusters)
    if (length(invalid_clusters) > 0) {
      warning(sprintf("The following clusters in cluster_order were not found in %s: %s",
                      group.by, paste(invalid_clusters, collapse = ", ")))
    }
  }
  
  obj_name <- deparse(substitute(seurat_object))
  
  # Prepare save path and create directory if needed
  final_save_path <- NULL
  if (!is.null(save_path)) {
    final_save_path <- save_path
    save_dir_path <- dirname(final_save_path)
    if (!dir.exists(save_dir_path)) {
      dir.create(save_dir_path, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("Created directory: %s", save_dir_path))
    }
  } else if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      message(sprintf("Created directory: %s", save_dir))
    }
    file_name <- paste0(obj_name, "_stack_vln.png")
    final_save_path <- file.path(save_dir, file_name)
  }
  
  features <- unique(features)
  available_features <- base::intersect(features, rownames(seurat_object[[assay]]))
  
  if(length(available_features) == 0) {
    stop(sprintf("None of the requested features were found in assay '%s'. Requested: %s",
                 assay, paste(features, collapse=", ")))
  }
  
  if(length(available_features) < length(features)){
    missing <- base::setdiff(features, available_features)
    warning(sprintf("The following %d gene(s) were not found in assay '%s' and will be skipped: %s", 
                    length(missing), assay, paste(missing, collapse=", ")))
    features <- available_features 
  }
  
  df_raw <- FetchData(seurat_object, vars = c(group.by, features), layer = "data")
  
  df_long <- df_raw %>%
    pivot_longer(cols = all_of(features), names_to = "gene", values_to = "expression") %>%
    dplyr::rename(cluster = !!sym(group.by)) %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(!is.na(cluster))
  
  # 統計量
  df_stat <- df_long %>%
    group_by(cluster, gene) %>%
    summarise(
      mean_expr = mean(expression, na.rm = TRUE),   
      median_expr = median(expression, na.rm = TRUE), 
      .groups = "drop"
    )
  
  # Build cluster tree using Seurat's BuildClusterTree
  dendrogram_method <- match.arg(dendrogram_method)
  
  if (dendrogram_method == "features") {
    # Use only the features specified for plotting
    seurat_object <- BuildClusterTree(
      seurat_object,
      assay = assay,
      features = features,  # Use validated features
      slot = "data",
      reorder = FALSE,
      verbose = FALSE
    )
  } else if (dendrogram_method == "dims") {
    # Use dimensionality reduction space
    if (is.null(ndim)) {
      dims_for_tree <- 1:30
      message("ndim not specified, using default dims = 1:30 for dendrogram calculation")
    } else {
      # Validate ndim is a positive integer
      if (!is.numeric(ndim) || length(ndim) != 1 || ndim <= 0 || ndim != as.integer(ndim)) {
        stop("ndim must be a single positive integer")
      }
      dims_for_tree <- 1:ndim
    }
    
    # Validate that the reduction exists
    if (!reduction_for_tree %in% names(seurat_object@reductions)) {
      stop(sprintf("Reduction '%s' not found in Seurat object. Available reductions: %s",
                   reduction_for_tree, paste(names(seurat_object@reductions), collapse = ", ")))
    }
    
    seurat_object <- BuildClusterTree(
      seurat_object,
      assay = assay,
      dims = dims_for_tree,
      reduction = reduction_for_tree,
      reorder = FALSE,
      verbose = FALSE
    )
  } else {
    # Use all variable features
    seurat_object <- BuildClusterTree(
      seurat_object,
      assay = assay,
      slot = "data",
      reorder = FALSE,
      verbose = FALSE
    )
  }
  
  # Extract the dendrogram
  hc <- Tool(seurat_object, slot = "BuildClusterTree")
  
  # Continue with existing dendro_data plotting code
  dendro_data <- ggdendro::dendro_data(hc)
  p_dendro <- ggplot(ggdendro::segment(dendro_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() + 
    theme(plot.margin = margin(b = -10, unit = "pt")) 
  
  if (is.null(cluster_order)) {
    final_cluster_order <- as.character(hc$labels[hc$order])
  } else {
    final_cluster_order <- as.character(unique(cluster_order))
  }
  
  # 順序付け
  df_plot <- df_long %>%
    left_join(df_stat %>% select(cluster, gene, median_expr), by = c("cluster", "gene"))
  
  df_plot$cluster <- factor(df_plot$cluster, levels = final_cluster_order)
  df_plot$gene <- factor(df_plot$gene, levels = features) 
  
  if(any(is.na(df_plot$cluster))) warning("Some clusters became NA after reordering.")
  
  # 描画
  p_violin <- ggplot(df_plot, aes(x = cluster, y = expression, fill = median_expr)) +
    geom_violin(scale = "width", width = 1.0, trim = TRUE, linewidth = 0.1, color = "black") + 
    facet_grid(rows = vars(gene), scales = "free_y", switch = "y") +
    scale_fill_gradient(low = color_low, high = color_high, name = "Median\nExpr") +
    theme_classic() +
    theme(
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 10, face = "bold"),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black")
    )
  
  p_final <- p_dendro / p_violin + plot_layout(heights = plot_heights)
  
  if (!is.null(final_save_path)) {
    calc_height <- max(5, length(features) * 0.4 + 2)
    tryCatch({
      ggsave(
        filename = final_save_path,
        plot = p_final,
        width = plot_width,
        height = calc_height,
        dpi = 300,
        bg = "white",
        limitsize = FALSE
      )
      message(sprintf("Successfully saved plot to: %s (Width: %d in, Height: %.1f in)", 
                      final_save_path, plot_width, calc_height))
    }, error = function(e) {
      warning(sprintf("Failed to save plot to '%s': %s", final_save_path, e$message))
    })
  }
  
  return(p_final)
}