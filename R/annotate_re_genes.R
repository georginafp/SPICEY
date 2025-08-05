#' Extract promoter regions annotated gene symbols from a TxDb and AnnotationDbi object
#' @inheritParams annotate_with_coaccessibility
#' @return A GRanges object with the chromosomes, start and end positions
#' of defined specie promoter regions together with the official gene symbol
#' stored in the `gene_id` metadata column.
get_promoters <- function(txdb,
                          annot_dbi,
                          upstream,
                          downstream,
                          protein_coding_only = TRUE) {
  proms <- GenomicFeatures::promoters(GenomicFeatures::genes(txdb),
    upstream = upstream,
    downstream = downstream
  )
  entrez_ids <- names(proms)
  gene_info <- AnnotationDbi::select(
    annot_dbi,
    keys = entrez_ids,
    columns = c("SYMBOL", "GENETYPE"),
    keytype = "ENTREZID"
  ) |>
    filter(!duplicated(ENTREZID))

  if (protein_coding_only) {
    gene_info <- gene_info[gene_info$GENETYPE == "protein-coding", ]
  }

  keep_ids <- gene_info$ENTREZID
  proms <- proms[entrez_ids %in% keep_ids]
  mcols(proms)$gene_id <- gene_info$SYMBOL[match(names(proms), gene_info$ENTREZID)]

  return(proms)
}





#' Annotates regulatory elements (e.g., ATAC-seq peaks) to the nearest gene
#' based on distance to the transcription start site (TSS), using a \code{TxDb}
#' reference and optional gene annotations from \code{org.*.db} packages.
#' @inheritParams annotate_with_coaccessibility
#' @return A \code{data.frame} of peaks annotated to its nearest gene, with columns:
#'   \itemize{
#'     \item \code{distanceToTSS}: distance to the nearest TSS
#'     \item \code{gene_id}: Official gene symbol of the nearest gene (e.g., GAPDH)
#'     \item \code{annotation}: \code{"Promoter"} or \code{"Distal"} based on distance}
#' @export
#' @examples
#' library(dplyr)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' 
#' data(atac)
#' peaks <- unique(unlist(atac)[,c("region_id")])
#' 
#' annotation_near <- annotate_with_nearest(
#'   peaks = peaks,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db,
#'   protein_coding_only = TRUE,
#'   verbose = TRUE,
#'   add_tss_annotation = FALSE,
#'   upstream = 2000,
#'   downstream = 2000
#' )
annotate_with_nearest <- function(peaks,
                                  txdb,
                                  annot_dbi,
                                  protein_coding_only = TRUE,
                                  verbose = TRUE,
                                  add_tss_annotation = FALSE,
                                  upstream,
                                  downstream) {
  if (verbose) message("Annotating regulatory elements to nearest gene...")
  if (inherits(peaks, "data.frame") && !inherits(peaks, "GRanges")) {
    peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  }
  alt_chroms <- grep("_alt|random|fix|Un", seqlevels(peaks), value = TRUE)
  peaks <- GenomeInfoDb::keepSeqlevels(peaks, setdiff(seqlevels(peaks), alt_chroms), pruning.mode = "coarse")

  ref_anno <- extract_gene_peak_annotations(
      peaks, txdb, annot_dbi,
      protein_coding_only,
      upstream, downstream, verbose)

  ref_anno <- ref_anno |>
    data.frame(row.names = NULL) |>
    tidyr::separate(peak,
      into = c("chr", "start", "end"),
      sep = "-", convert = TRUE
    ) |>
    dplyr::distinct() |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  nearest_hits <- GenomicRanges::distanceToNearest(peaks, ref_anno)
  dist <- rep(NA_integer_, length(peaks))
  genes <- rep(NA_character_, length(peaks))
  dist[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance
  genes[queryHits(nearest_hits)] <- ref_anno$gene_id[subjectHits(nearest_hits)]
  peaks$distanceToTSS <- dist
  peaks$gene_id <- genes
  peaks$annotation <- ifelse(!is.na(dist) & abs(dist) <= 2000, "Promoter", "Distal")

  if (add_tss_annotation) {
    peaks <- annotate_tss(peaks,
                          txdb = txdb,
                          annot_dbi = annot_dbi,
                          protein_coding_only = protein_coding_only)
  }

  return(as.data.frame(peaks))
}




#' Annotate regulatory elements overlapping transcription start sites (TSS)
#'
#' Identifies regulatory elements that overlap precisely defined transcription
#' start sites (TSS) and assigns the corresponding gene symbols to these elements.
#' @inheritParams annotate_with_coaccessibility
#' @return A \code{data.frame} containing the input regulatory elements with added columns:
#' \describe{
#'   \item{\code{in_TSS}}{Logical indicating whether the element overlaps a TSS.}
#'   \item{\code{TSS_gene}}{Gene symbol of the overlapping TSS, or \code{NA} if none.}}
#' @importFrom GenomicRanges makeGRangesFromDataFrame promoters findOverlaps mcols
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors queryHits subjectHits
annotate_tss <- function(txdb, peaks, annot_dbi, protein_coding_only) {
  proms <- get_promoters(
    txdb = txdb,
    annot_dbi = annot_dbi,
    upstream = 0,
    downstream = 1,
    protein_coding_only = protein_coding_only
  )
  hits <- findOverlaps(peaks, proms)
  GenomicRanges::mcols(peaks)$in_TSS <- FALSE
  GenomicRanges::mcols(peaks)$TSS_gene <- NA_character_
  GenomicRanges::mcols(peaks)$in_TSS[queryHits(hits)] <- TRUE
  GenomicRanges::mcols(peaks)$TSS_gene[queryHits(hits)] <- GenomicRanges::mcols(proms)$gene_id[subjectHits(hits)]
  peaks <- peaks |> data.frame()
  return(peaks)
}





#' Overlap peaks with gene promoters to obtain gene annotations
#'
#' Identifies overlaps between a set of peaks and promoter regions,
#' optionally restricted to protein-coding genes.
#' @inheritParams annotate_with_coaccessibility
#' @return A \code{data.frame} with:
#' \describe{
#'   \item{gene_id}{Identifier of the gene. This must be official gene symbols (e.g., GAPDH)}
#'   \item{peak}{Unique identifier of the region (e.g., chr1-5000-5800}}
extract_gene_peak_annotations <- function(peaks,
                                          txdb,
                                          annot_dbi,
                                          protein_coding_only = TRUE,
                                          upstream,
                                          downstream,
                                          verbose = FALSE) {
  proms <- get_promoters(
    txdb = txdb,
    annot_dbi = annot_dbi,
    upstream = upstream,
    downstream = downstream,
    protein_coding_only = protein_coding_only
  )

  ols <- GenomicRanges::findOverlaps(peaks, proms)
  gene_peak_df <- tibble::tibble(
    gene_id = mcols(proms)$gene_id[subjectHits(ols)],
    peak = granges_to_string(peaks[queryHits(ols)])
  ) |>
    dplyr::distinct() |>
    data.frame()

  return(gene_peak_df)
}




#' Annotate Cicero co-accessibility links with genes
#' Assigns gene annotations to Cicero co-accessibility links by
#' matching peaks to promoter-associated genes.
#' @inheritParams annotate_with_coaccessibility
#' @param gene_peak_anno A \code{data.frame} returned by
#' \code{extract_gene_peak_annotations()},
#' with columns \code{gene_id} and \code{peak}.
#' @return A \code{data.frame} with:
#' \describe{
#'   \item{peak}{Unique identifier of the distal peak (e.g., chr1-5000-5800}
#'   \item{promoter_peak}{Unique identifier of the promoter-associated peak (e.g., chr1-5000-5800}
#'   \item{coaccess}{Cicero co-accessibility score.}
#'   \item{gene_id}{Identifier of the gene. This must be official gene symbols (e.g., GAPDH)}}
annotate_links_with_genes <- function(links_df, gene_peak_anno) {
  joined_links <- links_df |>
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
      gene_id = gene_name1
    ),
    dplyr::transmute(
      dplyr::filter(joined_links, !is.na(gene_name2)),
      peak = Peak1,
      promoter_peak = Peak2,
      coaccess,
      gene_id = gene_name2
    )
  ) |>
    dplyr::distinct()
  return(annot_links)
}







#' Annotate peaks with co-accessible genes using Cicero links
#' Links peaks to genes based on Cicero co-accessibility with promoters or TSSs.
#' @param peaks A \code{GRanges} or \code{data.frame} of peaks with at least the following columns:
#' \describe{
#'   \item{seqnames}{Chromosome name of the regulatory region (e.g., \code{"chr1"}).}
#'   \item{start}{Start coordinate of the peak.}
#'   \item{end}{End coordinate of the peak.}
#'   \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800})}}
#' @param links_df A \code{data.frame} with Cicero links. Must contain columns:
#' \code{Peak1}, \code{Peak2}, and \code{coaccess}.
#' @param protein_coding_only Logical; restrict to protein-coding genes (default TRUE).
#' @param txdb \code{TxDb} object for genome annotation (required if annotation requested).
#' @param annot_dbi \code{AnnotationDbi} object for gene ID mapping (required if annotation requested).
#' @param add_tss_annotation Logical; annotate regulatory elements overlapping TSS (default FALSE).
#' If TRUE, use +/- 1bp TSS instead of +/-2kb promoter.
#' @param upstream Single integer value indicating the number of bases upstream
#' from the TSS (transcription start sites) (default 2000kb).
#' @param downstream Single integer values indicating the number of bases downstream
#' from the TSS (transcription start sites) (default 2000kb).
#' @param verbose Logical; print messages (default TRUE).
#' @return A \code{data.frame} with the original metadata columns from \code{peaks},
#' along with an added \code{gene_id} column containing the symbol of the co-accessible gene.
#' Peaks with no gene annotation will have \code{NA} in the \code{gene_id} field.
#' @export
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' 
#' data(atac)
#' data(cicero_links)
#' 
#' peaks <- unique(unlist(atac)[,c("region_id")])
#' annotation_coacc <- annotate_with_coaccessibility(
#'   peaks = peaks,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   links_df = cicero_links,
#'   annot_dbi = org.Hs.eg.db,
#'   protein_coding_only = TRUE,
#'   verbose = TRUE,
#'   add_tss_annotation = FALSE,
#'   upstream = 2000,
#'   downstream = 2000
#' )
annotate_with_coaccessibility <- function(peaks,
                                          txdb,
                                          links_df,
                                          annot_dbi,
                                          protein_coding_only = TRUE,
                                          verbose = TRUE,
                                          add_tss_annotation = FALSE,
                                          upstream,
                                          downstream) {
  if (verbose) {
    message("Annotating with co-accessibility")
  }
  if (inherits(peaks, "data.frame") && !inherits(peaks, "GRanges")) {
    peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks,
      keep.extra.columns = TRUE
    )
  }
  alt <- grep("_alt|random|fix|Un", seqlevels(peaks), value = TRUE)
  peaks <- GenomeInfoDb::keepSeqlevels(peaks, setdiff(seqlevels(peaks), alt),
    pruning.mode = "coarse"
  )
  ref_anno <- extract_gene_peak_annotations(
      peaks, txdb, annot_dbi,
      protein_coding_only,
      upstream, downstream, verbose
    )

  links_anno <- annotate_links_with_genes(links_df, ref_anno)

  if (add_tss_annotation) {
    peaks <- annotate_tss(peaks,
                          txdb = txdb,
                          annot_dbi = annot_dbi,
                          protein_coding_only = protein_coding_only)
  }

  result <- as.data.frame(peaks) |>
    dplyr::left_join(dplyr::select(links_anno, region_id = peak, gene_id),
      by = "region_id"
    ) |>
    dplyr::distinct()

  return(result)
}
