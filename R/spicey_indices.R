#' Compute cell type specificity scores from single-cell RNA and/or ATAC data
#'
#' Computes:
#' \itemize{
#'   \item GETSI (Gene Expression Tissue Specificity Index) from single-cell RNA-seq differential expression data.
#'   \item RETSI (Regulatory Element Tissue Specificity Index) from single-cell ATAC-seq differential accessibility data.}
#' Either RNA or ATAC input must be provided.
#' @inheritParams SPICEY
#' @param diff Single \code{data.frame} with differential results, with required columns:
#'   \describe{
#'     \item{id}{Identifier of the gene/region (e.g., gene symbol, Ensembl ID). 
#'          The name of this column should be provided in argument \code{id}}
#'     \item{avg_log2FC}{Average log2 fold-change for the gene in that cell type.}
#'     \item{p_val}{Raw p-value for the differential test.}
#'     \item{p_val_adj}{Adjusted p-value (e.g., FDR-corrected).}
#'     \item{cell_type}{Cell type or cluster label.}.
#'     }
#' @return A data frame with specificity scores for each feature:
#'   \describe{
#'     \item{score (GETSI or RETSI)}{Specificity score (weighted log2FC).}
#'     \item{norm_entropy}{Shannon entropy-based specificity measure.}
compute_spicey_index <- function(diff = NULL,
                                 id = NULL) {
  index <- specificity_index(diff, group_col = id)
  entropy <- entropy_index(index, group_col = id)
  
  final <- index |> 
    dplyr::left_join(entropy, by = id) |>
    dplyr::select(-c(avg_FC, max_FC, norm_FC, weight, entropy))
  
  return(final)
}

#' Calculate specificity scores for grouped features
#'
#' This function computes a specificity index for different features
#' (e.g., genes or regions) based on differential expression or accessibility data.
#' It rescales fold-change values and weights them by significance
#' to quantify how specific a feature's activity is to a particular cell type.
#' @param da A data.frame containing differential results with at least the following columns:
#'   \describe{
#'     \item{avg_log2FC}{Average log2 fold-change of the feature (gene or region).}
#'     \item{p_val}{Raw p-value from the differential test.}
#'     \item{cell_type}{Cell type or cluster label.}
#'     \item{\code{[group_col]}}{Column containing the feature identifier (e.g., gene_id or region)
#'   The **name of this column must match the value passed to the `group_col` argument**}}
#' @param group_col A string specifying the name of the column in \code{da} that
#' identifies each feature, such as \code{gene_id} for genes or \code{region} for ATAC peaks.
#' @return A data.frame identical to the input but with additional columns:
#'   \describe{
#'     \item{avg_FC}{Fold-change converted from log2 scale.}
#'     \item{p_val_adj}{FDR-adjusted p-values computed across all entries.}
#'     \item{max_FC}{Maximum fold-change observed within each feature group.}
#'     \item{weight}{Normalized significance weight derived from adjusted p-values.}
#'     \item{norm_FC}{Fold-change normalized by maximum fold-change in the group.}
#'     \item{score}{Specificity score computed as the product of normalized fold-change signficantly weighted.}}
#' @importFrom stats p.adjust
specificity_index <- function(da, group_col) {
  stopifnot(group_col %in% colnames(da))
  index <- da |>
    dplyr::mutate(
      avg_FC = 2^avg_log2FC,
      p_val_adj = p.adjust(p_val, method = "fdr")) |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::mutate(
      max_FC = max(avg_FC, na.rm = TRUE),
      p_val_adj = ifelse(p_val_adj == 0,
                         min(p_val_adj[p_val_adj > 0], na.rm = TRUE),
                         p_val_adj),
      weight = scales::rescale(-log10(p_val_adj), to = c(0, 1))) |>
    dplyr::group_by(cell_type, .data[[group_col]]) |>
    dplyr::mutate(
      norm_FC = avg_FC / max_FC,
      score = norm_FC * weight) |>
    dplyr::ungroup() |>
    data.frame()
  
  return(index)
}

#' Calculate normalized shannon-entropy of specificity scores
#'
#' Computes the entropy of specificity scores (RETSI or GETSI) across cell types
#' Entropy quantifies how evenly a feature's activity is distributed among cell types,
#' yielding scores from 0 to 1, where  where values close to 1 indicate
#' widespread distribution across cell types, and values near 0 denote dominating
#' distribution towards one cell type.
#' @param spec_df A data.frame containing the computed specificity scores
#' containing at least the following columns:
#'   \describe{
#'     \item{cell_type}{Cell type or cluster label.}
#'     \item{score}{Specificity score for each feature in each cell type.}
#'     \item{\code{[group_col]}}{Column containing the feautre identifier (e.g., gene_id or region)
#'   The **name of this column must match the value passed to the `group_col` argument**}}
#' @param group_col A string specifying the name of the column in \code{da} that
#' identifies each feature, such as \code{gene_id} for genes or \code{region} for ATAC peaks.
#' @return A data.frame with one row per feature, containing:
#'   \describe{
#'     \item{group_col}{Feature identifier.}
#'     \item{entropy}{Raw Shannon entropy computed from specificity scores.}
#'     \item{norm_entropy}{Normalized entropy score (1 - exp(-entropy))
#'     bounded between 0 and 1, where lower values indicate higher specificity.}}
entropy_index <- function(spec_df, group_col) {
  spec_df <- data.frame(spec_df)
  n_celltypes <- length(unique(spec_df$cell_type))
  spec_df <- spec_df |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::mutate(
      prob = score / sum(score, na.rm = TRUE),
      entropy_component = -prob * log2(prob)) |>
    dplyr::summarise(
      entropy = sum(entropy_component, na.rm = TRUE),
      norm_entropy = 1 - exp(-entropy),
      .groups = "drop") |>
    data.frame()
  
  return(spec_df)
}
