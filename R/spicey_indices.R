#' Compute specificity index for grouped features
#'
#' Computes a specificity score for features (e.g., genes or regions)
#' by rescaling fold-change and significance values.
#' @param df A data frame containing differential accessibility
#' or expression data.
#' @param group_col A string indicating the name of the column
#' to group by (e.g., "gene_id" or "region").
#' @importFrom dplyr mutate group_by ungroup
#' @importFrom stats p.adjust
#' @importFrom scales rescale
#' @return A data frame with an additional column `score`
#' representing the specificity index.
specificity_index <- function(df, group_col) {
  stopifnot(group_col %in% colnames(df))

  df |>
    dplyr::mutate(
      avg_FC = 2^avg_log2FC,
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::mutate(
      max_FC = max(avg_FC, na.rm = TRUE),
      p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0],
                                             na.rm = TRUE), p_val_adj),
      weight = rescale(-log10(p_val_adj), to = c(0, 1))
    ) |>
    dplyr::group_by(cell_type, .data[[group_col]]) |>
    dplyr::mutate(
      norm_FC = avg_FC / max_FC,
      score = norm_FC * weight
    ) |>
    dplyr::ungroup() |>
    data.frame()
}


#' Compute normalized entropy of specificity scores
#'
#' General function to compute entropy given a dataframe, grouping column, and score column.
#'
#' @param df Data.frame with columns for group_col, cell_type, and score_col.
#' @param group_col Character string name of grouping column ("gene_id" or "region").
#' @return Data.frame with group_col and norm_entropy columns.
#' @importFrom dplyr group_by mutate summarise
entropy_index <- function(df, group_col) {
  df <- data.frame(df)
  n_celltypes <- length(unique(df$cell_type))

  df |>
    dplyr::group_by(.data[[group_col]]) |>
    dplyr::mutate(
      prob = score / sum(score, na.rm = TRUE),
      entropy_component = -prob * log2(prob)
    ) |>
    dplyr::summarise(
      entropy = sum(entropy_component, na.rm = TRUE),
      norm_entropy = 1 - exp(-entropy),
      .groups = "drop"
    ) |>
    data.frame()
}



#' Compute GETSI and/or RETSI specificity scores
#'
#' This function computes GETSI scores for differential expression/accessibility
#' data on scRNA-seq and/or scATAC-seq respectively, based on inputs provided.
#' Both inputs are optional, but at least one must be given.
#' @param rna A named list of Dataframe objects, typically one per cell type,
#'  containing RNA differential expression data with columns such as
#'  `avg_log2FC`, `p_val`, `gene_id` etc.
#' @param gene_id A string indicating the name of the column in each data frame
#' that contains gene identifiers (e.g., `"gene_id"` or `"gene"`).
#' Required if `rna` is provided.
#' @param atac A named list where each element corresponds to a
#' cell type. Each must include `avg_log2FC`, `p_val`, and genomic coordinates.
#' @importFrom stats p.adjust
#' @importFrom dplyr rename
#' @return If only RNA or ATAC is provided, returns a single data.frame
#' with scores and entropy.
#' If both provided, returns a list with elements getsi and retsi.
#' @export
#' @examples
#' # Load example data (included with the package)
#' data(rna)
#' data(atac)
#'
#' # Compute both GETSI and RETSI
#' spicey_al <- compute_spicey_index(rna = rna, atac = atac, gene_id = "gene_id")
#'
#' # Compute GETSI only
#' spicey_getsi <- compute_spicey_index(rna = rna, gene_id = "gene_id")
#'
#' # Compute RETSI only
#' spicey_retsi <- compute_spicey_index(atac = atac)
compute_spicey_index <- function(rna = NULL,
                         atac = NULL,
                         gene_id = NULL) {
  if (is.null(rna) && is.null(atac)) {
    stop("You must provide at least one of 'rna' or 'atac'.")
  }

  result <- list()

  # GETSI
  if (!is.null(rna)) {
    if (is.null(gene_id)) stop("You must specify 'gene_id' for RNA input.")
    message("Computing GETSI...")
    rna_df <- combine_gex_da(rna, gene_id = gene_id)
    getsi <- specificity_index(rna_df, group_col = gene_id)

    message("Computing entropy on GETSI...")
    entropy_df <- entropy_index(getsi, group_col = gene_id)

    getsi_final <- getsi |>
      dplyr::left_join(entropy_df, by = gene_id) |>
      dplyr::select(-c(avg_FC, max_FC, norm_FC, weight, entropy)) |>
      dplyr::rename(GETSI = score,
                    GETSI_entropy = norm_entropy)

    result$getsi <- getsi_final
  }

  # RETSI
  if (!is.null(atac)) {
    message("Computing RETSI...")
    if (inherits(atac[[1]], "GRanges")) {
      atac_df <- unlist(GenomicRanges::GRangesList(atac)) |> data.frame()
    } else {
      atac_df <- do.call(rbind, atac) |> data.frame()
    }
    atac_df$region <- with(atac_df, paste0(seqnames, ":", start, "-", end))
    retsi <- specificity_index(atac_df, group_col = "region")

    message("Computing entropy on RETSI...")
    entropy_df <- entropy_index(retsi, group_col = "region")

    retsi_final <- retsi |>
      dplyr::left_join(entropy_df, by = "region") |>
      dplyr::select(-c(avg_FC, max_FC, norm_FC, weight, entropy)) |>
      dplyr::rename(RETSI = score,
                    RETSI_entropy = norm_entropy)

    result$retsi <- retsi_final
  }

  # Return single df or list depending on input
  if (!is.null(rna) && is.null(atac)) {
    return(result$getsi)
  } else if (is.null(rna) && !is.null(atac)) {
    return(result$retsi)
  } else {
    return(result)
  }
}
