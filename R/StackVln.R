#' Stacked Violin Plot (Scanpy Style)
#'
#' SeuratオブジェクトからScanpy風のスタックバイオリンプロットを描画します。
#' デンドログラムによるクラスタ順序の自動調整と、画像保存機能を備えています。
#'
#' @param seurat_object Seurat object containing the assay data.
#' @param features Vector of gene names to plot.
#' @param group.by Column name in meta.data for grouping (default: "seurat_clusters").
#' @param cluster_order Optional vector to specify cluster order manually.
#' @param assay Assay to use (default: "SCT").
#' @param color_low Color for low expression (default: "white").
#' @param color_high Color for high expression (default: "#BD2130").
#' @param plot_heights Vector (len=2) relative heights of dendrogram and violin (default: c(1, 9)).
#' @param plot_width Width of the saved plot in inches (default: 10).
#' @param save_path Full path to save the file. Overrides save_dir.
#' @param save_dir Directory to save the file using a default filename.
#'
#' @return A patchwork object combining the dendrogram and violin plots.
#' @export
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
    assay = "SCT",
    color_low = "white",
    color_high = "#BD2130", 
    plot_heights = c(1, 9),
    plot_width = 10,
    save_path = NULL,
    save_dir = NULL 
) {
  
  obj_name <- deparse(substitute(seurat_object))
  
  final_save_path <- NULL
  if (!is.null(save_path)) {
    final_save_path <- save_path
  } else if (!is.null(save_dir)) {
    file_name <- paste0(obj_name, "_stack_vln.png")
    final_save_path <- file.path(save_dir, file_name)
  }
  
  features <- unique(features)
  available_features <- base::intersect(features, rownames(seurat_object[[assay]]))
  
  if(length(available_features) < length(features)){
    missing <- base::setdiff(features, available_features)
    message("Warning: The following genes were not found and skipped: ", paste(missing, collapse=", "))
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
      mean_expr = mean(expression),   
      median_expr = median(expression), 
      .groups = "drop"
    )
  
  # デンドログラム
  mat_mean <- df_stat %>%
    select(cluster, gene, mean_expr) %>%
    pivot_wider(names_from = cluster, values_from = mean_expr) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  d <- dist(t(mat_mean))
  hc <- hclust(d, method = "complete")
  
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
    ggsave(
      filename = final_save_path,
      plot = p_final,
      width = plot_width,
      height = calc_height,
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    message(paste("Saved plot to:", final_save_path, "(Width:", plot_width, "in)"))
  }
  
  return(p_final)
}