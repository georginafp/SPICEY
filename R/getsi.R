#' Combine Differential Expression Results Across Cell Types
#'
#' This function takes a list of differential expression data frames (from scRNA-seq)
#' for multiple cell types and combines them into a single annotated data frame.
#' Each data frame should contain gene-level statistics (e.g., log2FC, p-value).
#' @param rna_da A named list of data frames. Each element corresponds to a cell type
#'   and contains differential expression results with genes as row names.
#' @return A single data frame combining all input differential expression results,
#'   with added columns for gene ids and cell types.
#' @examples
#' # Assuming rna_da is a named list of DESeq2 or Seurat DE results per cell type:
#' combined_df <- combine_gex_da(rna_da)
#' @export
combine_gex_da <- function(rna_da) {
  gr_list_annot <- lapply(names(rna_da), function(cell_type) {
    df <- rna_da[[cell_type]]
    df$gene_id <- rownames(df)

    if (!"cell_type" %in% colnames(df)) {
      df$cell_type <- cell_type
    }
    df  # last evaluated expression returned
  })

  names(gr_list_annot) <- names(rna_da)
  gr_list_annot <- Filter(Negate(is.null), gr_list_annot)

  combined_df <- do.call(rbind, gr_list_annot)
  rownames(combined_df) <- NULL
  return(combined_df)
}



#' Compute GETSI scores on RNA data
#'
#' @param rna_da Dataframe object with RNA differential expression data (avg_log2FC, p_val, cell_type, gene_id (gene_id, hgnc, ensembl_id...))
#' @return Dataframe with computed GETSI scores (GETSI column)
#' @import dplyr
#' @export
getsi <- function(rna_da) {
  rna_da_gr <- combine_gex_da(rna_da)

  getsi_df <- data.frame(rna_da_gr) |>
    mutate(
      avg_FC = 2^(avg_log2FC),
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    group_by(gene_id) |>
    mutate(
      max_FC = max(avg_FC, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    group_by(cell_type) |>
    mutate(
      p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm=TRUE), p_val_adj),
      weight = scales::rescale(-log10(p_val_adj), to=c(0, 1))
    ) |>
    group_by(cell_type, gene_id) |>
    mutate(
      norm_FC = avg_FC / max_FC,
      GETSI = (norm_FC * weight)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(gene_id, everything()) |>
    data.frame()
  return(getsi_df)
}



#' Compute normalized entropy of GETSI scores
#'
#' @param getsi_df Dataframe with GETSI scores and gene_id column
#' @return Data frame with gene_id gene_id and norm_entropy columns
#' @import dplyr
#' @export
entropy_getsi <- function(getsi_df) {
  df <- data.frame(getsi_df)
  n_celltypes <- length(unique(df$cell_type))

  df <- df |>
    group_by(gene_id) |>
    mutate(
      GETSI_prob = GETSI / sum(GETSI, na.rm = TRUE),
      entropy_component = -GETSI_prob * log2(GETSI_prob)
    ) |>
    summarise(
      entropy = sum(entropy_component, na.rm = TRUE),
      # norm_entropy = entropy / log2(n_celltypes),
      norm_entropy = 1 - exp(-entropy),
      .groups = "drop"
    )
  return(df)
}


#' Compute GETSI scores and entropy from RNA-seq data
#'
#' Computes gene expression tissue specificity (GETSI) scores and
#' entropy from RNA-seq data. The input should be a named list of
#' Dataframe objects, each containing differential expression data for
#' one cell type.
#'
#' @param rna_da A named list of Dataframe objects, typically one per cell type,
#'   containing RNA differential expression data with columns such as
#'   avg_log2FC, p_val, gene_id, etc.
#'
#' @return A Dataframe object with GETSI scores and entropy values in the
#'   metadata columns.
#'
#' @export
spicey_getsi <- function(rna_da) {
  message("→ Computing GETSI scores...")
  getsi_df <- getsi(rna_da)

  message("→ Computing entropy on GETSI scores...")
  entropy_df <- entropy_getsi(getsi_df)

  # Join entropy info with getsi results
  getsi_final <- as.data.frame(getsi_df) |>
    dplyr::left_join(entropy_df, by = "gene_id") |>
    dplyr::select(-c(avg_FC,
                     max_FC, norm_FC,
                     weight, entropy))
  return(getsi_final)
}
