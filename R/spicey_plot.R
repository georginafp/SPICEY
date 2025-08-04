#' Prepare data for SPICEY heatmap
#'
#' This function filters and processes a gene–cell-type matrix for heatmap visualization.
#' It computes z-scores for a specified score column `score_col` (e.g., `RETSI`, `GETSI`),
#' selects the top `n` genes per cell type based on their z-scores, and returns
#' a summary matrix suitable for plotting with \code{\link{plot_heatmap}}.
#' @param df A data frame containing at least the columns: \code{gene_id}, \code{cell_type},
#'   and one score column (e.g., \code{RETSI} or \code{GETSI}).
#' @param score_col Character. Name of the score column to z-score and rank (e.g., \code{"RETSI"}).
#' @param top_n Integer. Number of top-ranked genes to retain per cell type.
#' @return A data frame with columns: \code{gene_id}, \code{cell_type}, and \code{z_score}.
#'   Each row corresponds to a gene-cell-type pair, with averaged z-score values across duplicates.
#' @seealso \code{\link{plot_heatmap}}, \code{\link{spicey_heatmap}}
#' @import dplyr
prepare_heatmap_data <- function(df, score_col, top_n) {
  df_filtered <- df |>
    filter(!is.na(gene_id), !is.na(.data[[score_col]])) |>
    mutate(z_score = scale(.data[[score_col]])[, 1])

  top_genes <- df_filtered |>
    group_by(cell_type) |>
    arrange(desc(z_score)) |>
    distinct(cell_type, gene_id, .keep_all = TRUE) |>
    slice_head(n = top_n) |>
    ungroup() |>
    pull(gene_id)

  return(df_filtered |>
           filter(gene_id %in% top_genes) |>
           group_by(gene_id, cell_type) |>
           summarise(z_score = mean(z_score), .groups = "drop"))
}




#' Plot a z-scored gene-by-cell-type heatmap
#'
#' Generates a heatmap using ggplot2 to visualize expression or accessibility
#' z-scores for genes across different cell types. Genes are ordered by their
#' highest-scoring cell type, and then by maximum z-score within that group.
#' @param df_z A data frame with z-scored values. Must contain:
#'   \describe{
#'     \item{\code{gene_id}}{Identifier of the gene. This must be official gene symbols (e.g., GAPDH)}
#'     \item{\code{cell_type}}{Cell type or cluster label (e.g., Acinar)}
#'     \item{\code{z_score}}{Numeric. Z-scored values of the specificity score (e.g., RETSI_z, GETSI_z)}}
#' @param title_text Character. Title of the heatmap.
#' @param fill_label Character. Legend label for the color scale.
#' @return A \code{ggplot2} object representing the heatmap.
#' @importFrom textshape column_to_rownames
#' @importFrom tidyr pivot_wider
#' @import ggplot2
#' @import dplyr
#' @seealso \code{\link{prepare_heatmap_data}}, \code{\link{spicey_heatmap}}
plot_heatmap <- function(df_z, title_text, fill_label) {
  wide_mat <- df_z |>
    pivot_wider(names_from = cell_type, values_from = z_score, values_fill = 0) |>
    column_to_rownames("gene_id") |>
    as.matrix()

  gene_order_df <- data.frame(
    gene_id = rownames(wide_mat),
    max_cell_type = colnames(wide_mat)[apply(wide_mat, 1, which.max)],
    max_score = apply(wide_mat, 1, max)
  ) |>
    arrange(max_cell_type, desc(max_score))

  df_z$gene_id <- factor(df_z$gene_id, levels = gene_order_df$gene_id)

  return(ggplot(df_z, aes(x = cell_type, y = gene_id, fill = z_score)) +
    geom_tile(color = "grey80", width = 0.95) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = 0.5) +
    scale_fill_gradientn(
      colours = c("white", "#E18C80", "#B86357", "#723D36", "#4F2A25"),
      name = fill_label,
      na.value = "grey90"
    ) +
    theme_gray(base_size = 7) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 4)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(x = "Cell type", y = "Gene", title = title_text))
}



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
#'     \item{\code{gene_id}}{Identifier of the gene. This must be official gene symbols (e.g., GAPDH)}
#'     \item{\code{cell_type}}{Cell type or cluster label (e.g., Acinar)}
#'     \item{\code{RETSI}}{Numeric. RETSI specificity scores (optional unless used).}
#'     \item{\code{GETSI}}{Numeric. GETSI specificity scores (optional unless used).}}
#' @param top_n Integer. Number of top-ranked genes to include per cell type (default: 5).
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
#' library(dplyr)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' library(SPICEY)
#' data(rna)
#' data(atac)
#' data(cicero_links)
#' retsi <- SPICEY(atac = atac, region_id = "region_id")
#' getsi <- SPICEY(rna = rna, gene_id = "gene_id")
#' both <- SPICEY(
#'   rna = rna,
#'   gene_id = "gene_id",
#'   atac = atac,
#'   region_id = "region_id"
#' )
#' peaks <- SPICEY:::.parse_input_diff(atac)
#' peaks <- peaks %>%
#'   tidyr::separate(region_id,
#'     into = c("chr", "start", "end"), sep = "-",
#'     convert = TRUE, remove = FALSE
#'   ) %>%
#'   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
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
#' spicey_coacc <- SPICEY(
#'   rna = rna,
#'   gene_id = "gene_id",
#'   atac = atac,
#'   region_id = "region_id",
#'   annotation = annotation_coacc
#' )
#' spicey_heatmap(spicey_coacc$linked, spicey_measure = "RETSI")
#' spicey_heatmap(spicey_coacc$linked, spicey_measure = "GETSI")
#' spicey_heatmap(spicey_coacc$linked, spicey_measure = "SPICEY", combined_zscore = FALSE)
#' spicey_heatmap(spicey_coacc$linked, spicey_measure = "SPICEY", combined_zscore = TRUE)
#' @seealso \code{\link{SPICEY}}, \code{\link{prepare_heatmap_data}}, \code{\link{plot_heatmap}}
#' @export
#' @import dplyr
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
    return(plot_heatmap(df_z, spicey_measure, paste0(spicey_measure, "\nz-score")))
  }
  if (combined_zscore) {
    df_combined <- df |>
      filter(!is.na(gene_id), !is.na(RETSI), !is.na(GETSI)) |>
      mutate(
        GETSI_z = scale(GETSI)[, 1],
        RETSI_z = scale(RETSI)[, 1],
        combined_score = (RETSI_z + GETSI_z) / 2
      )
    top_genes <- df_combined |>
      group_by(cell_type) |>
      arrange(desc(combined_score)) |>
      distinct(cell_type, gene_id, .keep_all = TRUE) |>
      slice_head(n = top_n) |>
      ungroup() |>
      pull(gene_id)
    heatmap_df <- df_combined |>
      filter(gene_id %in% top_genes) |>
      group_by(gene_id, cell_type) |>
      summarise(z_score = mean(combined_score), .groups = "drop")
    return(plot_heatmap(heatmap_df, "SPICEY", "SPICEY\nz-score"))
  } else {
    df_z1 <- prepare_heatmap_data(df, "RETSI", top_n)
    df_z2 <- prepare_heatmap_data(df, "GETSI", top_n)
    return(plot_grid(plot_heatmap(df_z1, "RETSI", "RETSI\nz-score"),
                     plot_heatmap(df_z2, "GETSI", "GETSI\nz-score"), nrow = 1))
  }
}
