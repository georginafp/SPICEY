#' Compute RETSI scores on ATAC peaks
#'
#' @param gr GRanges object with differential accessibility info (avg_log2FC, p_val, cell_type columns)
#' @return GRanges with computed RETSI scores (RETSI column)
#' @importFrom regioneR toGRanges
#' @import dplyr
#' @export
retsi <- function(gr) {
  n_clusters <- length(unique(gr$cell_type))
  df <- data.frame(gr) |>
    mutate(
      region = paste0(seqnames, ":", start, "-", end),
      avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
      avg_FC = 2^(avg_FC),
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    group_by(region) |>
    mutate(
      max_log2FC = max(avg_FC, na.rm = TRUE),
      max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)
    ) |>
    ungroup() |>
    mutate(
      p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj == 0, 1e-300, p_val_adj),
      weight = -log10(p_val_adj)
    ) |>
    group_by(cell_type, region) |>
    mutate(
      norm_log2FC = (avg_FC / max_log2FC),
      RETSI = (norm_log2FC * weight) / (n_clusters - 1),
      RETSI = log1p(RETSI)
    ) |>
    ungroup() |>
    dplyr::select(-c(width, strand)) |>
    data.frame()
  return(regioneR::toGRanges(df))
}

#' Compute normalized entropy of RETSI scores
#'
#' @param gr GRanges with RETSI scores and region column
#' @return Data frame with region and norm_entropy columns
#' @import dplyr
#' @export
entropy_retsi <- function(gr) {
  df <- data.frame(gr) |>
    group_by(region) |>
    mutate(RETSI_sum = sum(RETSI, na.rm = TRUE),
           RETSI_prob = RETSI / RETSI_sum,
           RETSI_prob = RETSI_prob / sum(RETSI_prob, na.rm = TRUE),
           entropy_component = -RETSI_prob * log2(RETSI_prob)) |>
    summarise(entropy = sum(entropy_component, na.rm = TRUE)) |>
    # mutate(norm_entropy = 1 - exp(-entropy)) |>
    mutate(
      H_min = min(entropy, na.rm = TRUE),
      H_max = max(entropy, na.rm = TRUE),
      norm_entropy = (entropy - H_min) / (H_max - H_min)) |>
    dplyr::select(-c(H_min, H_max))
  return(df)
}



#' Compute RETSI scores and entropy from ATAC-seq data
#'
#' @param atac GRanges or path to RDS file with ATAC-seq data
#'
#' @return GRanges object with RETSI scores and entropy
#' @export
spicey_retsi <- function(atac) {
  if (is.character(atac)) atac <- readRDS(atac)

  message("→ Computing RETSI scores...")
  retsi <- retsi(atac)
  entropy <- entropy_retsi(retsi)

  mcols(retsi)$region <- paste0(seqnames(retsi), ":", start(retsi), "-", end(retsi))
  df <- data.frame(retsi) |>
    dplyr::left_join(entropy, by = "region") |>
    dplyr::select(-width, -strand)
  return(GenomicRanges::GRanges(df))
}

