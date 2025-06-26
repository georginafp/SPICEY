#' Annotate co-accessible links with CCAN membership and return GInteractions
#'
#' @param links A data.frame of co-accessibility links with columns Peak1, Peak2, coaccess.
#' @param coaccess_cutoff_override Numeric, coaccessibility cutoff for CCAN generation (default 0.25).
#' @param tolerance_digits Integer, rounding precision for cutoff (default 2).
#' @param filter_promoter_distal Logical, whether to keep only Promoter-Distal links (default TRUE).
#' @param txdb TxDb object for peak annotation (required if filter_promoter_distal=TRUE).
#'
#' @return GInteractions object annotated with coaccess, CCAN1, CCAN2 metadata columns.
#' @export
annotate_links_with_ccans <- function(links,
                                      coaccess_cutoff_override = 0.25,
                                      tolerance_digits = 2,
                                      filter_promoter_distal = TRUE,
                                      txdb = NULL) {
  # Generate CCAN assignments
  ccan <- generate_ccans(links,
                         coaccess_cutoff_override = coaccess_cutoff_override,
                         tolerance_digits = tolerance_digits)

  # Join CCAN annotations
  links <- links |>
    dplyr::left_join(ccan |> dplyr::rename(CCAN1 = CCAN), by = c("Peak1" = "Peak")) |>
    dplyr::left_join(ccan |> dplyr::rename(CCAN2 = CCAN), by = c("Peak2" = "Peak"))

  if (filter_promoter_distal) {
    if (is.null(txdb)) stop("TxDb object required for promoter-distal filtering")

    annotate_unique_peaks <- function(peak_vec) {
      df_peaks <- tibble::tibble(peak = unique(peak_vec)) |>
        tidyr::separate(peak,
                        into = c("seqnames", "start", "end"),
                        sep = "-",
                        convert = TRUE,
                        remove = FALSE) |>
        data.frame()

      gr <- regioneR::toGRanges(df_peaks |> dplyr::select(seqnames, start, end))
      anno <- ChIPseeker::annotatePeak(gr, TxDb = txdb, verbose = FALSE) |>
        data.frame()

      tibble::tibble(
        peak = df_peaks$peak,
        distanceToTSS = anno$distanceToTSS,
        annotation = ifelse(abs(anno$distanceToTSS) <= 2000, "Promoter", "Distal")
      )
    }

    peak1_anno <- annotate_unique_peaks(links$Peak1)
    peak2_anno <- annotate_unique_peaks(links$Peak2)

    links <- links |>
      dplyr::left_join(peak1_anno, by = c("Peak1" = "peak")) |>
      dplyr::rename(distanceToTSS1 = distanceToTSS, annotation1 = annotation) |>
      dplyr::left_join(peak2_anno, by = c("Peak2" = "peak")) |>
      dplyr::rename(distanceToTSS2 = distanceToTSS, annotation2 = annotation) |>
      dplyr::filter(xor(annotation1 == "Promoter", annotation2 == "Promoter"))
  }

  # Convert to GInteractions
  links_gi <- InteractionSet::GInteractions(
    Signac::StringToGRanges(links$Peak1),
    Signac::StringToGRanges(links$Peak2),
    coaccess = links$coaccess,
    CCAN1 = links$CCAN1,
    CCAN2 = links$CCAN2
  )

  return(links_gi)
}

#' Extract CCAN membership from links for a given GRanges object
#'
#' @param links GInteractions with CCAN1 and CCAN2 metadata columns
#' @param gr GRanges to query
#' @param name_column Character, metadata column in gr to use as identifier (default "region")
#' @param split Character, "ccan" returns names per CCAN, "name" returns CCANs per name
#'
#' @return Named list split by CCAN or name as requested
#' @export
get_ccan <- function(links, gr, name_column = "region", split = c("ccan", "name")) {
  split <- match.arg(split)

  hits1 <- findOverlaps(gr, links, use.region = "first")
  hits2 <- findOverlaps(gr, links, use.region = "second")

  ccan_df <- rbind(
    data.frame(ccan = mcols(links[subjectHits(hits1)])$CCAN1,
               name = mcols(gr[queryHits(hits1)])[[name_column]]),
    data.frame(ccan = mcols(links[subjectHits(hits2)])$CCAN2,
               name = mcols(gr[queryHits(hits2)])[[name_column]])
  ) |>
    dplyr::filter(!is.na(ccan), !is.na(name)) |>
    unique()

  if (split == "ccan") {
    split(ccan_df$name, ccan_df$ccan)
  } else {
    split(ccan_df$ccan, ccan_df$name)
  }
}




#' Link regulatory elements to genes based on shared CCAN membership
#'
#' @param links GInteractions with CCAN metadata
#' @param re GRanges of regulatory elements with 'region' column
#' @param proms GRanges of promoters with 'symbol' column
#' @param name_links Character suffix for metadata column name (default "HPAP")
#'
#' @return GRanges of regulatory elements annotated with linked genes
#' @export
get_targets_links <- function(links, re, proms, name_links = "HPAP") {
  res <- get_ccan(links, re, name_column = "region", split = "name")
  re$CCAN <- res[re$region]

  genes <- get_ccan(links, proms, name_column = "symbol", split = "ccan")

  mcols(re)[[paste0("genes_", name_links)]] <- lapply(re$CCAN, function(x) {
    if (!is.null(x)) unlist(genes[as.character(x)]) else NA
  })

  return(re)
}




#' Annotate regulatory elements with gene targets using co-accessibility links
#'
#' @param links GInteractions with CCAN metadata
#' @param re GRanges of regulatory elements with 'region' column
#' @param txdb TxDb object for promoter extraction
#' @param name_links Character suffix for gene link metadata column (default "HPAP")
#'
#' @return GRanges of regulatory elements with un-nested gene annotations
#' @export
annotate_with_coaccessibility <- function(links, re, txdb, name_links = "HPAP") {
  proms <- get_promoters_protein_coding(txdb)
  re_annotated <- get_targets_links(links, re, proms, name_links = name_links)
  re_unnested <- tidyr::unnest(as.data.frame(re_annotated),
                               cols = !!rlang::sym(paste0("genes_", name_links))) |>
    data.frame() |>
    select(-c(CCAN)) |>
    regioneR::toGRanges()
  return(return(re_unnested))
}


#' Annotate regions with nearest gene and promoter information
#'
#' This function uses ChIPseeker to annotate peaks with distance to TSS and
#' classifies them as Promoter (<2kb) or Distal. It then adds the nearest gene symbol
#' using promoters from Biomart-annotated genes.
#'
#' @param retsi A GRanges object containing regulatory elements to annotate.
#' @param genes A GRanges object with gene coordinates and external_gene_name (optional).
#'              Defaults to biomart_genes()$gr.
#' @return A GRanges object with metadata columns: distanceToTSS, annotation, nearestGeneSymbol.
#' @importFrom ChIPseeker annotatePeak
#' @importFrom plyranges mutate
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
annotate_with_nearest <- function(retsi, genes = biomart_genes()$gr) {
  # Ensure gene ranges have gene symbol in metadata
  if (is.null(genes$external_gene_name)) {
    stop("Gene ranges must have `external_gene_name` in metadata.")
  }

  # Precompute promoters (TSS regions)
  promoters_gr <- GenomicRanges::promoters(genes, upstream = 0, downstream = 1)

  # Find nearest TSS and compute distance
  nearest_hits <- GenomicRanges::distanceToNearest(retsi, promoters_gr)

  # Extract distances and gene symbols
  retsi$distanceToTSS <- rep(NA_integer_, length(retsi))
  retsi$nearestGeneSymbol <- rep(NA_character_, length(retsi))
  retsi$distanceToTSS[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance
  retsi$nearestGeneSymbol[queryHits(nearest_hits)] <- genes$external_gene_name[subjectHits(nearest_hits)]

  # Annotate as Promoter vs. Distal
  retsi$annotation <- ifelse(
    !is.na(retsi$distanceToTSS) & abs(retsi$distanceToTSS) <= 2000,
    "Promoter",
    "Distal"
  )

  # Remove names (optional)
  names(retsi) <- NULL

  return(retsi)
}
