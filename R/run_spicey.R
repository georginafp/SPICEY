#' Run the SPICEY pipeline
#'
#' This function runs the SPICEY pipeline for scoring and linking regulatory elements (REs) and genes
#' based on single-cell ATAC-seq and RNA-seq data. It supports two linking strategies: co-accessibility-based
#' and nearest-gene annotation.
#'
#' @param atac_path Path to the RDS file containing ATAC-seq data as a `GRanges` or `data.frame` convertible to `GRanges`.
#' @param rna_path Path to the RDS file containing RNA-seq data as a `GRanges` or `data.frame` convertible to `GRanges`.
#' @param links_path Optional path to the RDS file containing co-accessibility links as a `data.frame`. Required if `linking_method = "coaccessibility"`.
#' @param linking_method Linking method to use; either `"coaccessibility"` or `"nearest"`.
#'
#' @return A `data.frame` containing annotated REs with RETSI/GETSI scores, entropy, and linked gene information.
#' The returned table will include RETSI scores and entropy for each region, and GETSI scores and entropy for the
#' linked genes according to the chosen linking strategy.
#'
#' @details
#' - RETSI: Regulatory Element Tissue Specificity Index.
#' - GETSI: Gene Expression Tissue Specificity Index.
#' - When using `"coaccessibility"` mode, the function reads a precomputed link file and performs RETSI/GETSI scoring, linking, and merging.
#' - When using `"nearest"` mode, the function links each RE to its nearest gene and merges GETSI scores accordingly.
#'
#' @examples
#' \dontrun{
#' # Co-accessibility-based linking
#' run_spicey(
#'   atac_path = "data/FINAL_ATAC.rds",
#'   rna_path = "data/FINAL_RNA.rds",
#'   links_path = "data/COACC_LINKS.rds",
#'   linking_method = "coaccessibility"
#' )
#'
#' # Nearest-gene-based linking
#' run_spicey(
#'   atac_path = "data/FINAL_ATAC.rds",
#'   rna_path = "data/FINAL_RNA.rds",
#'   linking_method = "nearest"
#' )
#' }
#'
#' @export

run_spicey <- function(atac_path, rna_path, coaccess_threshold, links_path = NULL, linking_method = c("coaccessibility", "nearest")) {
  linking_method <- match.arg(linking_method)

  message("→ Step 1: Reading input data...")
  atac <- read_input_file(atac_path)
  rna <- read_input_file(rna_path)

  if (linking_method == "coaccessibility") {
    if (is.null(links_path)) {
      stop("You must provide a 'links_path' when using linking_method = 'coaccessibility'")
    }
    links_df <- read_input_file(links_path)
  }

  message("→ Step 2: Computing RETSI and RETSI entropy...")
  atac_scored <- compute_retsi(atac)
  atac_entropy <- compute_entropy_retsi(atac_scored)
  mcols(atac_scored)$region <- paste0(seqnames(atac_scored), ":", start(atac_scored), "-", end(atac_scored))
  atac_scored_df <- as.data.frame(atac_scored) |>
    left_join(atac_entropy, by = "region") |>
    select(-c("width", "strand")) |>
    toGRanges()

  message("→ Step 3: Computing GETSI and GETSI entropy...")
  rna_scored <- compute_getsi(rna)
  rna_entropy <- compute_entropy_getsi(rna_scored)
  rna_scored_df <- as.data.frame(rna_scored) |>
    left_join(rna_entropy, by = "symbol") |>
    select(-c("width", "strand")) |>
    toGRanges()

  if (linking_method == "coaccessibility") {
    message("→ Step 4: Linking REs and genes via co-accessibility...")
    links <- make_links(links_df, coaccess_threshold)

    gr_links <- annotate_with_coaccessibility(
      links = links,
      retsi = atac_scored_df,
      getsi = rna_scored_df,
      name_column_peaks = "region",
      name_column_genes = "symbol"
    )

    getsi_df <- as.data.frame(rna_scored_df) |>
      select(symbol, GETSI, cell_type, norm_entropy) |>
      rename(GETSI_entropy = norm_entropy)

    gr_links <- gr_links |>
      data.frame() |>
      rename(RETSI_entropy = norm_entropy) |>
      left_join(getsi_df, by = c("genes_coacc" = "symbol", "cell_type")) |>
      select(c(seqnames, start, end, cell_type, annotation, distanceToTSS,
               nearestGeneSymbol, genes_coacc, region, RETSI, RETSI_entropy, GETSI,
               GETSI_entropy))

  } else if (linking_method == "nearest") {
    message("→ Step 4: Linking REs and genes via nearest gene annotation...")
    gr_links <- annotate_with_nearest(
      retsi = atac_scored_df,
      getsi = rna_scored_df
    ) |> data.frame() |>
      rename(nearestGeneSymbol  = genes_nearest) |>
      select(c(seqnames, start, end, cell_type, annotation, distanceToTSS,
               nearestGeneSymbol, region, RETSI, RETSI_entropy, GETSI,
               GETSI_entropy))
  }

  final_df <- gr_links |>
    data.frame() |>
    dplyr::select(-dplyr::any_of(c("width", "strand"))) |>
    regioneR::toGRanges()

  message("🌶️️SPICEY pipeline completed successfully.")
  return(final_df)
}
