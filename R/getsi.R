#' Compute GETSI scores on RNA data
#'
#' @param gr GRanges object with RNA differential expression data (avg_log2FC, p_val, cell_type, symbol)
#' @return GRanges with computed GETSI scores (GETSI column)
#' @importFrom regioneR toGRanges
#' @import dplyr
#' @export
getsi <- function(gr) {
  n_clusters <- length(unique(gr$cell_type))
  df <- data.frame(gr) |>
    mutate(
      region = paste0(seqnames, ":", start, "-", end),
      avg_FC = ifelse(is.na(avg_log2FC), 1e-6, avg_log2FC),
      avg_FC = 2^(avg_FC),
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    group_by(symbol) |>
    mutate(
      max_log2FC = max(avg_FC, na.rm = TRUE),
      max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)
    ) |>
    ungroup() |>
    mutate(
      p_val_adj = ifelse(is.na(p_val_adj) | p_val_adj == 0, 1e-300, p_val_adj),
      weight = -log10(p_val_adj)
    ) |>
    group_by(cell_type, symbol) |>
    mutate(
      norm_log2FC = avg_FC / max_log2FC,
      GETSI = (norm_log2FC * weight) / (n_clusters - 1),
      GETSI = log1p(GETSI)
    ) |>
    ungroup() |>
    dplyr::select(-c(width, strand)) |>
    data.frame()
  return(regioneR::toGRanges(df))
}
#' Compute normalized entropy of GETSI scores
#'
#' @param gr GRanges with GETSI scores and symbol column
#' @return Data frame with symbol and norm_entropy columns
#' @import dplyr
#' @export
entropy_getsi <- function(gr) {
  df <- data.frame(gr) |>
    group_by(symbol) |>
    mutate(GETSI_sum = sum(GETSI, na.rm = TRUE),
           GETSI_prob = GETSI / GETSI_sum,
           GETSI_prob = GETSI_prob / sum(GETSI_prob, na.rm = TRUE),
           entropy_component = -GETSI_prob * log2(GETSI_prob)) |>
    summarise(entropy = sum(entropy_component, na.rm = TRUE)) |>
    # mutate(norm_entropy = 1 - exp(-entropy)) |>
    mutate(
      H_min = min(entropy, na.rm = TRUE),
      H_max = max(entropy, na.rm = TRUE),
      norm_entropy = (entropy - H_min) / (H_max - H_min)) |>
    dplyr::select(-c(H_min, H_max))
  return(df)
}


#' Compute GETSI scores and entropy from RNA-seq data
#'
#' @param rna GRanges or path to RDS file with RNA-seq data
#'
#' @return GRanges object with GETSI scores and entropy
#' @export
spicey_getsi <- function(rna) {
  if (is.character(rna)) rna <- readRDS(rna)

  message("→ Computing GETSI scores...")
  getsi <- getsi(rna)
  entropy <- entropy_getsi(getsi)

  df <- as.data.frame(getsi) |>
    dplyr::left_join(entropy, by = "symbol") |>
    dplyr::select(-width, -strand)
  return(GenomicRanges::GRanges(df))
}
