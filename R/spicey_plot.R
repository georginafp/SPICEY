#' SPICEY heatmap for gene specificity across cell types
#'
#' Visualizes gene-level specificity scores (RETSI and/or GETSI) across cell types
#' using a z-scored heatmap representation. Depending on the chosen mode, the function
#' can display either RETSI or GETSI scores independently, or compute and visualize
#' a combined SPICEY score (mean z-score of RETSI and GETSI).
#' If \code{spicey_measure = "SPICEY"} and \code{combined_zscore = TRUE}, RETSI and GETSI
#' scores are scaled, averaged, and shown in a unified heatmap. Otherwise, separate
#' heatmaps are produced for RETSI and GETSI, respectively.
#' @param df A data frame with at least the following columns:
#'   \describe{
#'     \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH}).
#'     If you only have ATAC data, link to nearest gene (\code{\link{annotate_with_nearest}})
#'     or using coaccessibility (\code{\link{annotate_with_coaccessibility}}).}
#'     \item{\code{cell_type}}{Cell type or cluster label (e.g., \code{Acinar})}
#'     \item{\code{RETSI}}{Numeric. RETSI specificity scores (optional unless used).}
#'     \item{\code{GETSI}}{Numeric. GETSI specificity scores (optional unless used).}}
#' @param top_n Integer. Number of top-ranked genes to include per cell type (default \code{"5"})
#' @param spicey_measure Character. Score type to visualize. Must be one of the following:
#'   \describe{
#'     \item{\code{"RETSI"}}{Only RETSI will be plotted.}
#'     \item{\code{"GETSI"}}{Only GETSI will be plotted.}
#'     \item{\code{"SPICEY"}}{Both RETSI and GETSI are used (requires both columns)}}
#' @param combined_zscore Logical. Only relevant if \code{spicey_measure = "SPICEY"}.
#'   If \code{TRUE}, a single heatmap of mean RETSI/GETSI z-score is generated.
#'   If \code{FALSE}, two heatmaps are produced side by side (RETSI and GETSI).
#' @return A \code{ggplot2} object, or a \code{patchwork} layout if two heatmaps are returned.
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#'
#' data(rna)
#' data(atac)
#' data(cicero_links)
#'
#' # Obtain annotatin with coaccessibility
#' peaks <- unique(unlist(atac)[, c("region_id")])
#' annotation_coacc <- annotate_with_coaccessibility(
#'   peaks = peaks,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   links_df = cicero_links,
#'   annot_dbi = org.Hs.eg.db,
#'   protein_coding_only = TRUE,
#'   verbose = TRUE,
#'   add_tss_annotation = FALSE,
#'   upstream = 2000,
#'   downstream = 2000
#' )
#'
#' # Obtain linked SPICEY measures
#' spicey_coacc <- SPICEY(
#'   rna = rna,
#'   atac = atac,
#'   annotation = annotation_coacc
#' )
#'
#' # Make plots
#' retsi <- spicey_coacc$RETSI |> dplyr::left_join(annotation_coacc, by = c("region_id"))
#' spicey_heatmap(retsi, spicey_measure = "RETSI")
#'
#' spicey_heatmap(spicey_coacc$GETSI, spicey_measure = "GETSI")
#'
#' spicey_heatmap(spicey_coacc$linked, spicey_measure = "SPICEY", combined_zscore = FALSE)
#'
#' spicey_heatmap(spicey_coacc$linked, spicey_measure = "SPICEY", combined_zscore = TRUE)
#' @seealso \code{\link{SPICEY}}, \code{\link{prepare_heatmap_data}}, \code{\link{plot_heatmap}}
#' @export
#' @importFrom cowplot plot_grid
spicey_heatmap <- function(df,
                           top_n = 5,
                           spicey_measure = c("RETSI", "GETSI", "SPICEY"),
                           combined_zscore = FALSE) {
  spicey_measure <- match.arg(spicey_measure)

  if (!all(c("gene_id", "cell_type") %in% colnames(df))) {
    stop("Input data must contain 'gene_id' and 'cell_type'.")
  }

  has_RETSI <- "RETSI" %in% colnames(df)
  has_GETSI <- "GETSI" %in% colnames(df)
  if (spicey_measure == "RETSI" && !has_RETSI) stop("Missing RETSI column.")
  if (spicey_measure == "GETSI" && !has_GETSI) stop("Missing GETSI column.")
  if (spicey_measure == "SPICEY" && (!has_RETSI || !has_GETSI)) {
    stop("SPICEY requires both RETSI and GETSI.")
  }
  if (combined_zscore && spicey_measure != "SPICEY") {
    stop("combined_zscore = TRUE is only valid with spicey_measure = 'SPICEY'")
  }

  if (spicey_measure %in% c("RETSI", "GETSI")) {
    df_z <- prepare_heatmap_data(df, spicey_measure, top_n)
    final_plot <- plot_heatmap(df_z, spicey_measure, paste0(spicey_measure, "\nz-score"))

  } else if (spicey_measure == "SPICEY" & combined_zscore) {
    df_combined <- df |>
      dplyr::filter(!is.na(gene_id), !is.na(RETSI), !is.na(GETSI)) |>
      dplyr::mutate(
        GETSI_z = scale(GETSI)[, 1],
        RETSI_z = scale(RETSI)[, 1],
        combined_score = (RETSI_z + GETSI_z) / 2)

    heatmap_df <- prepare_heatmap_data(df_combined, "combined_score", top_n)
    final_plot <- plot_heatmap(heatmap_df, "SPICEY", "SPICEY\nz-score")

  } else if (spicey_measure == "SPICEY" & !combined_zscore) {
    df_z1 <- prepare_heatmap_data(df, "RETSI", top_n)
    df_z2 <- prepare_heatmap_data(df, "GETSI", top_n)

    final_plot <- cowplot::plot_grid(
      plot_heatmap(df_z1, "RETSI", "RETSI\nz-score"),
      plot_heatmap(df_z2, "GETSI", "GETSI\nz-score"),
      nrow = 1, align = "h")
  }
  return(final_plot)
}

#' Prepare data for SPICEY heatmap
#'
#' Filters and processes a gene–cell-type matrix for heatmap visualization.
#' Computes z-scores for a specified score column, selects the top `n` genes per cell type,
#' and returns a summary matrix suitable for plotting.
#' @param df A data frame with at least: \code{gene_id}, \code{cell_type}, and a score column.
#' @param score_col Character. Name of the score column to z-score and rank.
#' @param top_n Integer. Number of top-ranked genes per cell type.
#' @return A data frame with: \code{gene_id}, \code{cell_type}, and \code{z_score}.
#' @seealso \code{\link{plot_heatmap}}, \code{\link{spicey_heatmap}}
prepare_heatmap_data <- function(df, score_col, top_n) {
  df_filtered <- df |>
    dplyr::filter(!is.na(gene_id), !is.na(.data[[score_col]])) |>
    dplyr::mutate(z_score = scale(.data[[score_col]])[, 1])

  top_genes <- df_filtered |>
    dplyr::group_by(cell_type) |>
    dplyr::arrange(dplyr::desc(z_score)) |>
    dplyr::distinct(cell_type, gene_id, .keep_all = TRUE) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup() |>
    dplyr::pull(gene_id)

  df_final <- df_filtered |>
    dplyr::filter(gene_id %in% top_genes) |>
    dplyr::group_by(gene_id, cell_type) |>
    dplyr::summarise(z_score = mean(z_score), .groups = "drop")
  return(df_final)
}


#' Plot a z-scored gene-by-cell-type heatmap
#'
#' Generates a heatmap using ggplot2 to visualize expression or accessibility
#' z-scores for genes across different cell types. Genes are ordered by their
#' highest-scoring cell type, and then by maximum z-score within that group.
#' @param df_z A data frame with z-scored values. Must contain:
#'   \describe{
#'     \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH})}.
#'     \item{\code{cell_type}}{Cell type or cluster label (e.g., \code{Acinar})}
#'     \item{\code{z_score}}{Numeric. Z-scored values of the specificity score (e.g., \code{RETSI_z}, \code{GETSI_z})}}
#' @param title_text Character. Title of the heatmap.
#' @param fill_label Character. Legend label for the color scale.
#' @return A \code{ggplot2} object representing the heatmap.
#' @seealso \code{\link{prepare_heatmap_data}}, \code{\link{spicey_heatmap}}
plot_heatmap <- function(df_z, title_text, fill_label) {
  wide_mat <- df_z |>
    tidyr::pivot_wider(names_from = cell_type, values_from = z_score, values_fill = 0) |>
    tibble::column_to_rownames(var = "gene_id") |>
    as.matrix()

  gene_order_df <- data.frame(
    gene_id = rownames(wide_mat),
    max_cell_type = colnames(wide_mat)[apply(wide_mat, 1, which.max)],
    max_score = apply(wide_mat, 1, max)
  ) |>
    dplyr::arrange(max_cell_type, dplyr::desc(max_score))

  df_z$gene_id <- factor(df_z$gene_id, levels = gene_order_df$gene_id)

  plot <- ggplot2::ggplot(df_z, ggplot2::aes(x = cell_type, y = gene_id, fill = z_score)) +
    ggplot2::geom_tile(color = "grey75", width = 0.95) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_fixed(ratio = 0.5) +
    ggplot2::scale_fill_gradientn(
      colours = c("white", "#d3bfc1", "#b29092", "#926164", "#664346"),
      name = fill_label,
      na.value = "grey85"
    ) +
    ggplot2::theme_gray(base_size = 7) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 4)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      plot.margin = grid::unit(c(1, 1, 1, 1), "mm") # top, right, bottom, left
    ) +
    ggplot2::labs(x = "Cell type", y = "Gene", title = title_text)

  return(plot)
}
