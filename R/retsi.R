#' Compute RETSI scores for differential chromatin accessibility data
#'
#' This function calculates the Regulatory Element Tissue Specificity Index (RETSI)
#' from a list of differential accessibility results (either as GRanges or data.frames).
#'
#' @param atac_da A named list where each element corresponds to a cell type and contains either:
#'   \itemize{
#'     \item A \code{GRanges} object with differential peaks and metadata columns: \code{avg_log2FC}, \code{p_val}
#'     \item A \code{data.frame} with the same structure, including a \code{cell_type} column
#'   }
#'
#' @return A data.frame with RETSI scores per region and cell type. Includes columns:
#'   \code{region}, \code{RETSI}, and supporting metadata.
#'
#' @import dplyr
#' @export
retsi <- function(atac_da) {
  # Input handling
  if (inherits(atac_da[[1]], "GRanges")) {
    atac_df <- unlist(GenomicRanges::GRangesList(atac_da)) |>
      data.frame()
  } else {
    atac_df <- do.call(rbind, atac_da) |> data.frame()
  }

  # RETSI computation
  atac_df |>
    dplyr::mutate(
      region = paste0(seqnames, ":", start, "-", end),
      avg_FC = 2^avg_log2FC,
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    dplyr::group_by(region) |>
    dplyr::mutate(max_FC = max(avg_FC, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::group_by(cell_type) |>
    dplyr::mutate(
      p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE), p_val_adj),
      weight = scales::rescale(-log10(p_val_adj), to = c(0, 1))
    ) |>
    dplyr::group_by(cell_type, region) |>
    dplyr::mutate(
      norm_FC = avg_FC / max_FC,
      RETSI = norm_FC * weight
    ) |>
    dplyr::ungroup() |>
    data.frame()
}



#' Compute entropy of RETSI scores per regulatory region
#'
#' Calculates normalized entropy based on RETSI scores across cell types,
#' providing a measure of cell-type specificity per region.
#'
#' @param retsi A data.frame output from \code{retsi()} containing \code{region},
#'   \code{cell_type}, and \code{RETSI} columns.
#'
#' @return A data.frame with columns: \code{region} and \code{norm_entropy}.
#'
#' @import dplyr
#' @export
entropy_retsi <- function(retsi) {
  df <- data.frame(retsi)
  n_celltypes <- length(unique(df$cell_type))

  df |>
    dplyr::group_by(region) |>
    dplyr::mutate(
      RETSI_prob = RETSI / sum(RETSI, na.rm = TRUE),
      entropy_component = -RETSI_prob * log2(RETSI_prob)
    ) |>
    dplyr::summarise(
      entropy = sum(entropy_component, na.rm = TRUE),
      norm_entropy = 1 - exp(-entropy),
      .groups = "drop"
    ) |>
    data.frame()
}


#' Compute RETSI scores and entropy for scATAC-seq regulatory elements
#'
#' This is a wrapper function that takes a list of per-cell-type differential accessibility results
#' (as GRanges or data.frames), computes RETSI scores, and adds a normalized entropy score
#' to quantify cell-type specificity of each regulatory region.
#'
#' @param atac_da A named list of GRanges or data.frames (one per cell type),
#'   each with differential accessibility metrics including \code{avg_log2FC} and \code{p_val}.
#'
#' @return A data.frame with RETSI scores and normalized entropy for each region,
#'   including \code{region}, \code{cell_type}, \code{RETSI}, and \code{norm_entropy}.
#' @export
spicey_retsi <- function(atac_da) {
  # Validate input
  if (!is.list(atac_da)) {
    stop("Input must be a list of GRanges or data.frames")
  }

  first_class <- class(atac_da[[1]])
  if (!first_class %in% c("GRanges", "data.frame")) {
    stop("Each element in the list must be either a GRanges or a data.frame")
  }

  message("→ Computing RETSI scores...")
  retsi_df <- retsi(atac_da)

  message("→ Computing entropy of RETSI scores...")
  entropy_df <- entropy_retsi(retsi_df)

  # Merge entropy and return cleaned result
  retsi_final <- retsi_df |>
    dplyr::left_join(entropy_df, by = "region") |>
    dplyr::select(-c(avg_FC,
                     max_FC, norm_FC,
                     weight, entropy))
  data.frame()
  return(retsi_final)
}
