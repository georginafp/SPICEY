#' Extract promoter regions with gene symbols
#'
#' Returns promoter regions from a TxDb object, annotated with gene symbols.
#' @param txdb A TxDb object.
#' @param annot_dbi An AnnotationDbi object (e.g., org.Hs.eg.db).
#' @param upstream,downstream Number of bases upstream/downstream of TSS (default: 2000).
#' @param protein_coding_only Logical. If TRUE, restrict to protein-coding genes.
#' @return A GRanges object with gene symbols in the `gene_id` metadata column.
#' @export
get_promoters <- function(txdb,
                          annot_dbi,
                          upstream,
                          downstream,
                          protein_coding_only = TRUE) {

  proms <- GenomicFeatures::promoters(GenomicFeatures::genes(txdb),
                                      upstream = upstream,
                                      downstream = downstream)
  entrez_ids <- names(proms)
  gene_info <- AnnotationDbi::select(
    annot_dbi,
    keys = entrez_ids,
    columns = c("SYMBOL", "GENETYPE"),
    keytype = "ENTREZID") |>
    filter(!duplicated(ENTREZID))

  if (protein_coding_only) {
    gene_info <- gene_info[gene_info$GENETYPE == "protein-coding", ]
  }

  keep_ids <- gene_info$ENTREZID
  proms <- proms[entrez_ids %in% keep_ids]
  mcols(proms)$gene_id <- gene_info$SYMBOL[match(names(proms), gene_info$ENTREZID)]

  return(proms)
}




#' Annotate regulatory elements with nearest gene (by TSS or promoter)
#'
#' @param peaks GRanges or data.frame of regulatory elements.
#' @param txdb TxDb object.
#' @param annot_dbi AnnotationDbi object (e.g., org.Hs.eg.db).
#' @param protein_coding_only Logical, whether to restrict to protein-coding genes.
#' @param keep_mito Logical, whether to retain mitochondrial genes.
#' @param verbose Logical, whether to display progress messages.
#' @param add_tss_annotation Logical, if TRUE, use precise TSS instead of broader promoter region.
#' @param upstream Integer, upstream window from TSS for promoter (only used if add_tss_annotation = FALSE).
#' @param downstream Integer, downstream window from TSS for promoter (only used if add_tss_annotation = FALSE).
#' @return A data.frame of input peaks annotated with nearest gene and distance to TSS.
annotate_with_nearest <- function(peaks,
                                  txdb,
                                  annot_dbi,
                                  protein_coding_only = TRUE,
                                  keep_mito = FALSE,
                                  verbose = TRUE,
                                  add_tss_annotation = FALSE,
                                  upstream,
                                  downstream) {

  if (verbose) message("Annotating regulatory elements to nearest gene...")

  # Coerce to GRanges if needed
  if (inherits(peaks, "data.frame") && !inherits(peaks, "GRanges")) {
    peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  }

  # Remove alternative contigs
  alt_chroms <- grep("_alt|random|fix|Un", seqlevels(peaks), value = TRUE)
  peaks <- GenomeInfoDb::keepSeqlevels(peaks, setdiff(seqlevels(peaks), alt_chroms), pruning.mode = "coarse")

  ref_anno <- if (add_tss_annotation) {
    extract_gene_peak_annotations(
      peaks, txdb, annot_dbi,
      protein_coding_only, keep_mito,
      upstream = 0, downstream = 1, verbose)
  } else {
    extract_gene_peak_annotations(
      peaks, txdb, annot_dbi,
      protein_coding_only, keep_mito,
      upstream, downstream,verbose)
  }

  ref_anno <- ref_anno |>
    tidyr::separate(peak, into = c("chr", "start", "end"),
                    sep = "-", convert = TRUE) |>
    dplyr::distinct() |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  # Find nearest gene
  nearest_hits <- GenomicRanges::distanceToNearest(peaks, ref_anno)
  dist <- rep(NA_integer_, length(peaks))
  genes <- rep(NA_character_, length(peaks))
  dist[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance
  genes[queryHits(nearest_hits)] <- ref_anno$gene_id[subjectHits(nearest_hits)]

  # Annotate
  peaks$distanceToTSS <- dist
  peaks$nearestGeneSymbol <- genes
  peaks$annotation <- ifelse(!is.na(dist) & abs(dist) <= 2000, "Promoter", "Distal")
  return(as.data.frame(peaks))
}







#' Extract gene-peak overlaps from promoter regions
#'
#' @param peaks A `GRanges` object with peaks (with region_id in metadata).
#' @param txdb A transcript database (TxDb object).
#' @param annot_dbi Annotation DB (e.g. `org.Hs.eg.db`).
#' @param protein_coding_only Logical. Whether to restrict to protein-coding genes.
#' @param keep_mito Logical. Whether to retain mitochondrial genes.
#' @param upstream Integer. Bases upstream of TSS to include in promoter.
#' @param downstream Integer. Bases downstream of TSS.
#' @param verbose Logical. Print messages.
#' @return A `data.frame` with `gene_id` and `peak` columns.
#' @keywords internal
extract_gene_peak_annotations <- function(peaks,
                                          txdb,
                                          annot_dbi,
                                          protein_coding_only = TRUE,
                                          keep_mito = FALSE,
                                          upstream,
                                          downstream,
                                          verbose = FALSE) {

  proms <- get_promoters(txdb = txdb,
                         annot_dbi = annot_dbi,
                         upstream = upstream,
                         downstream = downstream,
                         protein_coding_only = protein_coding_only)

  # Overlap peaks with promoters
  ols <- GenomicRanges::findOverlaps(peaks, proms)
  gene_peak_df <- tibble::tibble(
    gene_id = mcols(proms)$gene_id[subjectHits(ols)],
    peak = granges_to_string(peaks[queryHits(ols)])) |>
    dplyr::distinct() |>
    data.frame()

  return(gene_peak_df)
}





#' Annotate Cicero co-accessibility links with gene names
#'
#' @param links A data.frame with Cicero-style links (`Peak1`, `Peak2`, `coaccess`).
#' @param gene_peak_anno A data.frame with columns `gene_id` and `peak`.
#' @return A long-format data.frame with `peak`, `promoter_peak`, `coaccess`, and `gene_id`.
#' @keywords internal
annotate_links_with_genes <- function(links, gene_peak_anno) {
  joined_links <- links |>
    dplyr::left_join(gene_peak_anno, by = c("Peak1" = "peak")) |>
    dplyr::rename(gene_name1 = gene_id) |>
    dplyr::left_join(gene_peak_anno, by = c("Peak2" = "peak")) |>
    dplyr::rename(gene_name2 = gene_id)

  annot_links <- dplyr::bind_rows(
    dplyr::transmute(
      dplyr::filter(joined_links, !is.na(gene_name1)),
      peak = Peak2,
      promoter_peak = Peak1,
      coaccess,
      gene_id = gene_name1),

    dplyr::transmute(
      dplyr::filter(joined_links, !is.na(gene_name2)),
      peak = Peak1,
      promoter_peak = Peak2,
      coaccess,
      gene_id = gene_name2)) |>
    dplyr::distinct()

  return(annot_links)
}







#' Annotate peaks with co-accessible genes using Cicero links
#'
#' Links peaks to genes based on Cicero co-accessibility with promoters or TSSs.
#'
#' @param peaks A `GRanges` or `data.frame` of peaks with `region_id`.
#' @param txdb A `TxDb` object with transcript annotations.
#' @param links_df A `data.frame` with Cicero links (`Peak1`, `Peak2`, `coaccess`).
#' @param annot_dbi Annotation DB (e.g. `org.Hs.eg.db`).
#' @param protein_coding_only Logical. Keep only protein-coding genes.
#' @param keep_mito Logical. Keep mitochondrial genes.
#' @param verbose Logical. Verbose output.
#' @param coaccess_cutoff_override Numeric. Minimum coaccessibility (default = 0.25).
#' @param add_tss_annotation Logical. If TRUE, use ±1bp TSS instead of ±2kb promoter.
#' @param upstream,downstream Integers for promoter definition (default = 2000).
#' @return A `data.frame` of peaks annotated with `gene_id`.
#' @export
annotate_with_coaccessibility <- function(peaks,
                                          txdb,
                                          links_df,
                                          annot_dbi,
                                          protein_coding_only = TRUE,
                                          keep_mito = FALSE,
                                          verbose = TRUE,
                                          coaccess_cutoff_override = 0.25,
                                          add_tss_annotation = FALSE,
                                          upstream,
                                          downstream) {

  if (verbose) {
    message("Annotating with co-accessibility (cutoff: ", coaccess_cutoff_override, ")...")
  }

  if (inherits(peaks, "data.frame") && !inherits(peaks, "GRanges")) {
    peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks,
                                                     keep.extra.columns = TRUE)
  }

  alt <- grep("_alt|random|fix|Un", seqlevels(peaks), value = TRUE)
  peaks <- GenomeInfoDb::keepSeqlevels(peaks, setdiff(seqlevels(peaks), alt),
                                       pruning.mode = "coarse")

  ref_anno <- if (add_tss_annotation) {
    extract_gene_peak_annotations(
      peaks, txdb, annot_dbi,
      protein_coding_only, keep_mito,
      upstream = 0, downstream = 1, verbose = FALSE)
  } else {
    extract_gene_peak_annotations(
      peaks, txdb, annot_dbi,
      protein_coding_only, keep_mito,
      upstream, downstream, verbose)
  }

  links <- dplyr::filter(links_df, coaccess > coaccess_cutoff_override)
  links_anno <- annotate_links_with_genes(links, ref_anno)

  result <- as.data.frame(peaks) |>
    dplyr::left_join(dplyr::select(links_anno, region_id = peak, gene_id),
                     by = "region_id") |>
    dplyr::distinct()

  return(result)
}


