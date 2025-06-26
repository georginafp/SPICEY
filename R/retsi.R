#' Compute RETSI scores on ATAC peaks
#'
#' @param atac_da GRanges object with differential accessibility info (avg_log2FC, p_val, cell_type columns)
#' @return GRanges with computed RETSI scores (RETSI column)
#' @importFrom regioneR toGRanges
#' @import dplyr
#' @export
retsi <- function(atac_da) {
  atac_da_gr <- unlist(GRangesList(atac_da))
  retsi_gr <- atac_da_gr |>
    data.frame() |>
    dplyr::select(c(seqnames, start, end, everything())) |>
    dplyr::select(-c(width, strand)) |>
    regioneR::toGRanges()

  # Fix here: use atac_da_gr or retsi_gr as appropriate
  retsi <- data.frame(atac_da_gr) |>  # <-- corrected
    dplyr::mutate(
      region = paste0(seqnames, ":", start, "-", end),
      avg_FC = 2^(avg_log2FC),
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    dplyr::group_by(region) |>
    dplyr::mutate(
      max_log2FC = max(avg_FC, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(cell_type) |>
    dplyr::mutate(
      p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm=TRUE), p_val_adj),
      weight = scales::rescale(-log10(p_val_adj), to=c(0, 1))
    ) |>
    dplyr::group_by(cell_type, region) |>
    dplyr::mutate(
      norm_log2FC = avg_FC / max_log2FC,
      RETSI = (norm_log2FC * weight)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-c(width, strand)) |>
    data.frame() |>
    regioneR::toGRanges()

  return(retsi)
}


#' Compute normalized entropy of RETSI scores
#'
#' @param gr GRanges with RETSI scores and region column
#' @return Data frame with region and norm_entropy columns
#' @import dplyr
#' @export
entropy_retsi <- function(retsi) {
  df <- data.frame(retsi)
  n_celltypes <- length(unique(df$cell_type))

  entropy <- df |>
    group_by(region) |>
    mutate(
      RETSI_prob = RETSI / sum(RETSI, na.rm = TRUE),
      entropy_component = -RETSI_prob * log2(RETSI_prob)
    ) |>
    summarise(
      entropy = sum(entropy_component, na.rm = TRUE),
      # norm_entropy = entropy / log2(n_celltypes),
      norm_entropy = 1 - exp(-entropy),
      .groups = "drop"
    )
  return(entropy)
}



#' Compute RETSI scores from a List of GRanges representing scATAC-seq peaks per cell type
#'
#' @param atac_list List of GRanges objects, each containing peaks (regulatory elements) per cell type
#'
#' @return A GRanges object with RETSI scores and entropy annotations
#'
#' @details
#' The function expects a list of GRanges (e.g. one GRanges per cell type or condition).
#' It combines the list into a single GRanges, computes RETSI scores and entropy,
#' then returns a GRanges with added metadata columns.
#'
#' @examples
#' \dontrun{
#' atac_list <- list(celltype1 = gr1, celltype2 = gr2)
#' retsi_gr <- spicey_retsi(atac_list)
#' }
#' @export
spicey_retsi <- function(atac_da) {
  # Validate input is a list
  if (!is.list(atac_da)) {
    stop("Input must be a list of GRanges objects")
  }

  # Check all elements are GRanges
  if (!all(sapply(atac_da, inherits, "GRanges"))) {
    stop("All elements of the input list must be GRanges objects")
  }

  message("→ Computing RETSI scores...")
  retsi_gr <- retsi(atac_da)

  message("→ Computing entropy of RETSI scores...")
  entropy <- entropy_retsi(retsi_gr)

  # Add 'region' ID for joining
  mcols(retsi_gr)$region <- paste0(seqnames(retsi_gr), ":",
                                   start(retsi_gr), "-",
                                   end(retsi_gr))

  # Join entropy data.frame by region
  retsi_final <- data.frame(retsi_gr) |>
    dplyr::left_join(entropy, by = "region") |>
    dplyr::select(-c(width, strand, avg_FC,
                     max_log2FC, norm_log2FC,
                     weight, entropy)) |>
    regioneR::toGRanges()

  return(retsi_final)
}
