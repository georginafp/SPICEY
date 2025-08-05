#' Link RETSI regions to GETSI scores using gene-based association methods
#'
#' This function connects regulatory regions scored with RETSI
#' @param retsi A data.frame containing RETSI scores for chromatin accessibility regions,
#'   as returned by \code{compute_spicey_index()} using single-cell ATAC-seq differential accessibility data.
#'   Must include at least the following columns:
#'   \describe{
#'     \item{region_id}{Unique identifier of the region (e.g., chr1-5000-5800)}
#'     \item{cell_type}{Cell type or cluster label.}
#'     \item{RETSI}{RETSI value: cell-type specificity score}
#'     \item{norm_entropy}{Normalized Shannon entropy of RETSI}
#'   }
#' @param getsi A data.frame containing GETSI scores for genes,
#'   as returned by \code{compute_spicey_index()} using single-cell RNA-seq differential expression data.
#'   Must include at least the following columns:
#'   \describe{
#'     \item{gene_id}{Identifier of the gene. This must be official gene symbols (e.g., GAPDH)}
#'     \item{cell_type}{Cell type or cluster label.}
#'     \item{GETSI}{GETSI value: cell-type specificity score}
#'     \item{norm_entropy}{Normalized Shannon entropy of GETSI}
#'   }
#' @inheritParams SPICEY
#' @return A \code{data.frame} where each row represents a regulatory element–gene pair
#'   linked within a given cell type. The output includes:
#'   \describe{
#'     \item{region_id}{Unique identifier of the region (e.g., chr1-5000-5800)}
#'     \item{gene_id}{Identifier of the gene. This must be official gene symbols (e.g., GAPDH)}
#'     \item{cell_type}{Cell type or cluster in which the association is observed.}
#'     \item{RETSI}{RETSI score: regulatory element specificity in this cell type.}
#'     \item{RETSI_entropy}{Normalized shannon-entropy of RETSI (lower = more specific).}
#'     \item{GETSI}{GETSI score: gene expression specificity in this cell type.}
#'     \item{GETSI_entropy}{Normalized shannon-entropy of GETSI (lower = more specific).}
#'     \item{...}{Any additional columns from the original \code{retsi} and \code{getsi} inputs, suffixed
#'     with \code{_ATAC} and \code{_RNA} respectively (e.g., \code{avg_log2FC_ATAC}, \code{p_val_RNA}).}
#'   }
link_spicey <- function(retsi = NULL,
                        region_id = NULL,
                        getsi = NULL,
                        gene_id = NULL,
                        annotation = NULL) {
  keep_cols <- c("region_id", "cell_type", "gene_id", "distanceToTSS",
                 "annotation", "TSS_gene", "in_TSS")

  links <- retsi |>
    dplyr::inner_join(
      annotation |>
        dplyr::select(dplyr::any_of(keep_cols)),
      by = c(region_id)
    ) |>
    dplyr::inner_join(
      getsi,
      by = c(gene_id),
      suffix = c("_ATAC", "_RNA")
    )

  return(links)
}

