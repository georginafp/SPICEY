#' Annotates regulatory elements (e.g., ATAC-seq peaks) to the nearest gene
#'
#' based on distance to the transcription start site (TSS), using a \code{TxDb}
#' reference and gene annotations from \code{org.*.db} packages.
#' @inheritParams annotate_with_coaccessibility
#' @return A \code{data.frame} of peaks annotated to its nearest gene, with columns:
#' \item{distanceToTSS}{Distance to the nearest TSS}
#' \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH})}
#' \item{annotation}{\code{"Promoter"} or \code{"Distal"} based on distance}
#' @examples
#' library(dplyr)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#'
#' data(atac)
#' peaks <- unique(unlist(atac)[, c("region_id")])
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
#' @export
annotate_with_nearest <- function(peaks,
                                  txdb,
                                  annot_dbi,
                                  protein_coding_only = TRUE,
                                  verbose = TRUE,
                                  add_tss_annotation = FALSE,
                                  upstream,
                                  downstream) {
  if (verbose) message("Annotating regulatory elements to nearest gene...")
  peaks <- .standardize_peaks(peaks)

  ref_anno <- extract_gene_peak_annotations(
    peaks, txdb, annot_dbi, protein_coding_only, upstream, downstream, verbose)

  nearest_hits <- GenomicRanges::distanceToNearest(peaks, ref_anno)
  annotation <- peaks
  annotation$distanceToTSS <- NA_real_
  annotation$distanceToTSS[S4Vectors::queryHits(nearest_hits)] <-
    GenomicRanges::mcols(nearest_hits)$distance

  annotation$gene_id <- NA_character_
  annotation$gene_id[S4Vectors::queryHits(nearest_hits)] <-
    ref_anno$gene_id[S4Vectors::subjectHits(nearest_hits)]

  annotation$annotation <- ifelse(annotation$region_id %in% unique(ref_anno$region_id),
                                  "Promoter", "Distal")

  annotation <- GenomicRanges::mcols(
    if (add_tss_annotation) {
      .add_tss_annotation(annotation = annotation, peaks = peaks,
                          txdb = txdb, annot_dbi = annot_dbi,
                          protein_coding_only = protein_coding_only,
                          verbose = verbose)
    } else annotation
  ) |> data.frame()

  return(annotation)
}

#' Annotate peaks with co-accessible genes using Cicero links
#'
#' Links peaks to genes based on Cicero co-accessibility with promoters or TSSs.
#' @param peaks A \code{GRanges} or \code{data.frame} of peaks with at least the following columns:
#' \describe{
#'   \item{seqnames}{Chromosome name of the regulatory region (e.g., \code{"chr1"}). Only for data.frames.}
#'   \item{start}{Start coordinate of the peak. Only for data.frames.}
#'   \item{end}{End coordinate of the peak. Only for data.frames.}
#'   \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800})}
#'   }
#' @param links_df A \code{data.frame} with Cicero links. Must contain columns:
#' \code{Peak1}, \code{Peak2}, and \code{coaccess}.
#' @param protein_coding_only Logical; restrict to protein-coding genes (default TRUE).
#' @param txdb \code{TxDb} object for genome annotation (required if annotation requested).
#' @param annot_dbi \code{AnnotationDbi} object for gene ID mapping (required if annotation requested).
#' @param add_tss_annotation Logical; annotate regulatory elements overlapping TSS (default FALSE).
#' If TRUE, use +/- 1bp TSS.
#' @param upstream Single integer value indicating the number of bases upstream
#' from the TSS (transcription start sites) (default \code{"2000kb"}).
#' @param downstream Single integer values indicating the number of bases downstream
#' from the TSS (transcription start sites) (default \code{"2000kb"}).
#' @param verbose Logical; print messages (default TRUE).
#' @return A \code{data.frame} with the original metadata columns from \code{peaks},
#' along with an added \code{gene_id} column containing the symbol of the co-accessible gene.
#' Peaks with no gene annotation will have \code{NA} in the \code{gene_id} field.
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#'
#' data(atac)
#' data(cicero_links)
#'
#' peaks <- unique(unlist(atac)[, c("region_id")])
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
#' @export
annotate_with_coaccessibility <- function(peaks,
                                          txdb,
                                          links_df,
                                          annot_dbi,
                                          protein_coding_only = TRUE,
                                          verbose = TRUE,
                                          add_tss_annotation = FALSE,
                                          upstream,
                                          downstream) {
  if (verbose) message("Annotating with co-accessibility")
  peaks <- .standardize_peaks(peaks)

  ref_anno <- extract_gene_peak_annotations(
    peaks, txdb, annot_dbi, protein_coding_only, upstream, downstream, verbose) |>
    GenomicRanges::mcols() |>
    data.frame()

  joined_links <- links_df |>
    dplyr::left_join(ref_anno, by = c("Peak1" = "region_id")) |>
    dplyr::rename(gene_name1 = gene_id) |>
    dplyr::left_join(ref_anno, by = c("Peak2" = "region_id")) |>
    dplyr::rename(gene_name2 = gene_id)

  links_anno <- dplyr::bind_rows(
    dplyr::transmute(dplyr::filter(joined_links, !is.na(gene_name1)),
                     region_id = Peak2, promoter_peak = Peak1, coaccess, gene_id = gene_name1),
    dplyr::transmute(dplyr::filter(joined_links, !is.na(gene_name2)),
                     region_id = Peak1, promoter_peak = Peak2, coaccess, gene_id = gene_name2)) |>
    dplyr::distinct()

  annotation <- links_anno |>
    dplyr::select(region_id, gene_id)

  if (add_tss_annotation) {
    annotation <- .add_tss_annotation(annotation, peaks, txdb, annot_dbi,
                                      protein_coding_only, verbose)
  }
  return(annotation)
}

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
                                      downstream = downstream)
  entrez_ids <- names(proms)
  gene_info <- AnnotationDbi::select(
    annot_dbi, keys = entrez_ids, columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID") |>
    dplyr::distinct(ENTREZID, .keep_all = TRUE)

  if (protein_coding_only) gene_info <- dplyr::filter(gene_info, GENETYPE == "protein-coding")

  keep_ids <- gene_info$ENTREZID
  proms <- proms[entrez_ids %in% keep_ids]
  GenomicRanges::mcols(proms)$gene_id <- gene_info$SYMBOL[match(names(proms), gene_info$ENTREZID)]
  return(proms)
}

#' Standardizes peak input, ensuring that input peaks is a \code{GRanges} object,
#' removes alternate scaffolds (\code{_alt}, \code{random}, \code{fix}, \code{Un}),
#' and assigns region IDs as names.
#' @param peaks A \code{data.frame} with genomic coordinates and a \code{region_id}
#'   column, or a \code{GRanges} object with a \code{region_id} column.
#' @return A \code{GRanges} object with:
#' \describe{
#'   \item{seqnames, start, end}{Coordinates of the input regions.}
#'   \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800}).}
#'   \item{metadata columns}{Any additional columns present in the input.}
#' }
.standardize_peaks <- function(peaks) {
  if (inherits(peaks, "data.frame") && !inherits(peaks, "GRanges")) {
    peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  }
  alt <- grep("_alt|random|fix|Un", GenomeInfoDb::seqlevels(peaks), value = TRUE)
  peaks <- GenomeInfoDb::keepSeqlevels(peaks, setdiff(GenomeInfoDb::seqlevels(peaks), alt), pruning.mode = "coarse")
  names(peaks) <- peaks$region_id
  return(peaks)
}

#' Overlap peaks with gene promoters to obtain gene annotations
#'
#' Identifies overlaps between a set of peaks and promoter regions,
#' optionally restricted to protein-coding genes.
#' @inheritParams annotate_with_coaccessibility
#' @return A \code{GRanges with} with:
#' \describe{
#'   \item{seqnames, start, end}{Coordinates of the peak that overlaps with a gene promoter.}
#'   \item{region_id}{Unique identifier of the region (e.g., \code{chr1-5000-5800}).}
#'   \item{gene_id}{Identifier of the gene. Must be an official gene symbol (e.g., \code{GAPDH})}
#'   }
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
    protein_coding_only = protein_coding_only)

  ols <- GenomicRanges::findOverlaps(peaks, proms)

  promoter_peaks <- peaks[S4Vectors::queryHits(ols)]
  promoter_peaks$gene_id <- proms$gene_id[S4Vectors::subjectHits(ols)]
  return(promoter_peaks)
}


#' Add TSS annotation to peaks
#'
#' Identifies transcription start sites (TSS) overlapping with input peaks and
#' adds logical and gene ID columns to an annotation table.
#' @inheritParams annotate_with_coaccessibility
#' @param annotation A \code{data.frame} or \code{GRanges} containing peak annotations,
#'   including a \code{region_id} column.
#' @param peaks A \code{GRanges} object containing the original peak set.
#' @return A \code{data.frame} or \code{GRanges} (depending on input) with added columns:
#' \describe{
#'   \item{in_TSS}{Logical, \code{TRUE} if the peak overlaps a TSS.}
#'   \item{TSS_gene}{Gene symbol of the overlapping TSS, if any.}
#' }
.add_tss_annotation <- function(annotation, peaks, txdb, annot_dbi, protein_coding_only, verbose) {
  tss <- extract_gene_peak_annotations(
    peaks        = peaks,
    txdb         = txdb,
    annot_dbi    = annot_dbi,
    protein_coding_only = protein_coding_only,
    upstream     = 0,
    downstream   = 1,
    verbose      = verbose)

  annotation <- annotation |>
    dplyr::mutate(in_TSS = region_id %in% tss$region_id) |>
    dplyr::left_join(
      tss |> data.frame() |> dplyr::select(region_id, TSS_gene = gene_id),
      by = "region_id")
  return(annotation)
}
