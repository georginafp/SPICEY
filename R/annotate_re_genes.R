#' Annotate regulatory elements overlapping transcription start sites (TSS)
#'
#' This function identifies regulatory elements that directly overlap annotated
#' transcription start sites (TSS, defined as 1 bp ranges) and assigns the corresponding gene symbols.
#'
#' @param txdb A \code{TxDb} object used to extract gene coordinates.
#' @param re A \code{GRanges} object or data.frame of regulatory elements to annotate.
#' @param annot_dbi An \code{AnnotationDbi} object (e.g., \code{org.Hs.eg.db}) for gene symbol mapping.
#'
#' @return A \code{data.frame} of regulatory elements annotated with:
#' \itemize{
#'   \item \code{in_TSS}: Logical indicating whether the element overlaps a TSS.
#'   \item \code{TSS_gene}: Gene symbol of the overlapping TSS, if any.
#' }
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame promoters findOverlaps mcols
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors queryHits subjectHits
annotate_tss <- function(txdb, re, annot_dbi) {

  prom_tss <- promoters(genes(txdb), upstream=0, downstream=1)  # TSS as 1 bp ranges
  gene_ids <- GenomicRanges::mcols(prom_tss)$gene_id
  symbols <- mapIds(annot_dbi,
                    keys=gene_ids,
                    column="SYMBOL",
                    keytype="ENTREZID",
                    multiVals="first")
  GenomicRanges::mcols(prom_tss)$symbol <- symbols

  # Convert to GRanges if needed
  if (inherits(re, "data.frame") && !inherits(re, "GRanges")) {
    re <- GenomicRanges::makeGRangesFromDataFrame(re,
                                                  keep.extra.columns = TRUE)
  }

  hits <- findOverlaps(re, prom_tss)

  # Initialize columns
  GenomicRanges::mcols(re)$in_TSS <- FALSE
  GenomicRanges::mcols(re)$TSS_gene <- NA_character_

  # Mark TRUE and add gene symbol for those overlapping
  GenomicRanges::mcols(re)$in_TSS[queryHits(hits)] <- TRUE
  GenomicRanges::mcols(re)$TSS_gene[queryHits(hits)] <- GenomicRanges::mcols(prom_tss)$symbol[subjectHits(hits)]
  re <- re |> data.frame()
  return(re)
}




#' Annotate co-accessible links with CCAN membership and return GInteractions
#' @param links A data.frame of co-accessibility links with columns Peak1, Peak2, coaccess.
#' @param coaccess_cutoff_override Numeric, coaccessibility cutoff for CCAN generation (default 0.25).
#' @param tolerance_digits Integer, rounding precision for cutoff (default 2).
#' @param filter_promoter_distal Logical, whether to keep only Promoter-Distal links (default TRUE).
#' @param txdb TxDb object for peak annotation (required if filter_promoter_distal=TRUE).
#' @return GInteractions object annotated with coaccess, CCAN1, CCAN2 metadata columns.
annotate_links_with_ccans <- function(links,
                                      coaccess_cutoff_override = 0.25,
                                      tolerance_digits = 2,
                                      filter_promoter_distal = TRUE,
                                      txdb = NULL) {
  # Generate CCAN assignments
  ccan <- cicero::generate_ccans(links,
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
#' @param links GInteractions with CCAN1 and CCAN2 metadata columns
#' @param gr GRanges to query
#' @param name_column Character, metadata column in gr to use as identifier (default "region")
#' @param split Character, "ccan" returns names per CCAN, "name" returns CCANs per name
#' @return Named list split by CCAN or name as requested
#' @importFrom GenomicRanges findOverlaps mcols
#' @importFrom S4Vectors queryHits subjectHits
get_ccan <- function(links, gr, name_column = "region", split = c("ccan", "name")) {
  split <- match.arg(split)

  hits1 <- findOverlaps(gr, links, use.region = "first")
  hits2 <- findOverlaps(gr, links, use.region = "second")

  ccan_df <- rbind(
    data.frame(ccan = GenomicRanges::mcols(links[subjectHits(hits1)])$CCAN1,
               name = GenomicRanges::mcols(gr[queryHits(hits1)])[[name_column]]),
    data.frame(ccan = GenomicRanges::mcols(links[subjectHits(hits2)])$CCAN2,
               name = GenomicRanges::mcols(gr[queryHits(hits2)])[[name_column]])) |>
    dplyr::filter(!is.na(ccan), !is.na(name)) |>
    unique()

  if (split == "ccan") {
    split(ccan_df$name, ccan_df$ccan)
  } else {
    split(ccan_df$ccan, ccan_df$name)
  }
}




#' Link regulatory elements to genes based on shared CCAN membership
#' @param links GInteractions with CCAN metadata
#' @param re GRanges of regulatory elements with 'region' column
#' @param proms GRanges of promoters with 'gene_id' column
#' @return GRanges of regulatory elements annotated with linked genes
get_targets_links <- function(links, re, proms) {
  res <- get_ccan(links, re, name_column = "region", split = "name")
  re$CCAN <- res[re$region]
  genes <- get_ccan(links, proms, name_column = "gene_id", split = "ccan")
  GenomicRanges::mcols(re)$gene_coacc <- lapply(re$CCAN, function(x) {
    if (!is.null(x)) unlist(genes[as.character(x)]) else NA
  })
  return(re)
}




#' Annotate Regulatory Elements with Gene Targets Using Co-Accessibility Links
#'
#' This function annotates regulatory elements (REs) with their putative gene targets by integrating
#' co-accessibility links that include CCAN (co-accessibility networks) metadata. It extracts promoter
#' regions from a TxDb object, optionally filters for protein-coding genes, and links REs to genes
#' based on chromatin co-accessibility. Optionally, REs overlapping a transcription start site (TSS)
#' can also be annotated directly.
#'
#' @param links A \code{GInteractions} object representing co-accessibility links with CCAN metadata.
#' @param retsi A \code{GRanges} or \code{data.frame} object representing regulatory elements to annotate.
#'              If a data.frame, must contain columns sufficient to construct \code{GRanges} (e.g., seqnames, start, end).
#' @param txdb A \code{TxDb} object used to extract promoter regions and TSS annotations.
#' @param annot_dbi An \code{AnnotationDbi} annotation database object (e.g., \code{org.Hs.eg.db}) for mapping gene IDs to gene symbols and types.
#' @param protein_coding_only Logical indicating whether to restrict promoter regions to protein-coding genes only (default: \code{TRUE}).
#' @param keep_mito Logical indicating whether to keep mitochondrial (chrM/MT) chromosomes in promoter extraction (default: \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, informative messages are printed during processing (default: \code{TRUE}).
#' @param coaccess_cutoff_override Numeric; co-accessibility score cutoff to override default filtering in \code{annotate_links_with_ccans} (default: \code{0.25}).
#' @param filter_promoter_distal Logical; whether to filter links to retain only promoter-distal interactions (default: \code{TRUE}).
#' @param add_tss_annotation Logical; whether to annotate regulatory elements overlapping TSS directly (default: \code{TRUE}).
#'
#' @return A \code{data.frame} with regulatory elements annotated with:
#' \itemize{
#'   \item Distance to nearest TSS and promoter/distal classification
#'   \item Co-accessible gene targets
#'   \item Optional: TSS overlap and TSS gene symbol
#' }
#'
#' @importFrom GenomicRanges promoters distanceToNearest makeGRangesFromDataFrame mcols
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#' @importFrom S4Vectors queryHits subjectHits
#' @export
#' @examples
#' data(cicero_links)
#' data(atac)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' retsi <- compute_spicey_index(atac=atac)
#' # Annotate REs using co-accessibility links to predict gene targets
#' retsi_gene_coacc <- annotate_with_coaccessibility(
#'   links = cicero_links,
#'   retsi = retsi,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
#'   coaccess_cutoff_override = 0.25
#' )
annotate_with_coaccessibility <- function(links,
                                          retsi,
                                          txdb,
                                          annot_dbi,
                                          protein_coding_only = TRUE,
                                          keep_mito = FALSE,
                                          verbose = TRUE,
                                          coaccess_cutoff_override = 0.25,
                                          filter_promoter_distal = TRUE,
                                          add_tss_annotation = TRUE) {
  if (verbose) {
    message("Annotating regulatory elements to co-accessible genes...")
    message("Coaccessibility cutoff used: ", coaccess_cutoff_override)
  }

  # Convert to GRanges if needed
  if (inherits(retsi, "data.frame") && !inherits(retsi, "GRanges")) {
    retsi <- GenomicRanges::makeGRangesFromDataFrame(retsi,
                                                     keep.extra.columns = TRUE)
  }

  # Remove ALT contigs
  alt_chroms <- grep("_alt|random|fix|Un", seqlevels(retsi), value = TRUE)
  retsi <- GenomeInfoDb::keepSeqlevels(retsi,
                                       setdiff(seqlevels(retsi), alt_chroms),
                                       pruning.mode = "coarse")

  # Extract promoter regions
  proms <- get_promoters(txdb = txdb,
                         annot_dbi = annot_dbi,
                         keep_mito = keep_mito,
                         protein_coding_only = protein_coding_only,
                         verbose = verbose)

  # Compute distance to nearest TSS
  nearest_hits <- GenomicRanges::distanceToNearest(retsi, proms)
  retsi$distanceToTSS <- NA_integer_
  retsi$distanceToTSS[queryHits(nearest_hits)] <- GenomicRanges::mcols(nearest_hits)$distance

  # Annotate as Promoter vs. Distal
  retsi$annotation <- ifelse(
    !is.na(retsi$distanceToTSS) & abs(retsi$distanceToTSS) <= 2000,
    "Promoter",
    "Distal"
  )

  # Annotate links with CCANs
  annotated_links <- annotate_links_with_ccans(
    links = links,
    coaccess_cutoff_override = coaccess_cutoff_override,
    filter_promoter_distal = filter_promoter_distal,
    txdb = txdb
  )

  # Get target genes from co-accessibility
  re_annotated <- get_targets_links(
    links = annotated_links,
    re = retsi,
    proms = proms
  )

  # Unnest gene_coacc column and drop CCAN
  re_unnested <- tidyr::unnest(as.data.frame(re_annotated),
                               cols = "gene_coacc") |>
    dplyr::select(-c(CCAN)) |>
    data.frame()

  # Optionally add TSS overlap annotation
  if (add_tss_annotation) {
    re_unnested <- annotate_tss(re_unnested,
                                txdb = txdb,
                                annot_dbi = annot_dbi)
  }

  return(re_unnested)
}







#' Annotate regions with nearest gene and promoter information
#'
#' This function uses ChIPseeker-like logic to annotate peaks with distance to the nearest
#' TSS (transcription start site), classifies them as Promoter (<2kb) or Distal, and adds
#' the nearest gene ID using Biomart-annotated genes. Optionally, it can annotate whether a
#' region overlaps a TSS directly using 1-bp TSS coordinates.
#'
#' @param retsi A \code{GRanges} object or \code{data.frame} containing regulatory elements to annotate.
#' @param txdb A \code{TxDb} object used to extract promoter regions.
#' @param annot_dbi An \code{AnnotationDbi} annotation database object (e.g., \code{org.Hs.eg.db}) for mapping gene IDs to gene symbols and types.
#' @param protein_coding_only Logical indicating whether to restrict promoter regions to protein-coding genes only (default: \code{TRUE}).
#' @param keep_mito Logical indicating whether to keep mitochondrial (chrM/MT) chromosomes in promoter extraction (default: \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, informative messages are printed during processing (default: \code{TRUE}).
#' @param add_tss_annotation Logical; whether to also annotate whether a region overlaps a TSS (default: \code{TRUE}).
#'
#' @return A \code{data.frame} with the following columns:
#' \itemize{
#'   \item \code{distanceToTSS}: Distance to the nearest TSS
#'   \item \code{annotation}: Promoter/Distal classification based on distance
#'   \item \code{nearestGeneSymbol}: Symbol of the nearest gene
#'   \item \code{in_TSS} (optional): Logical indicating if RE overlaps a TSS
#'   \item \code{TSS_gene} (optional): Gene symbol for overlapping TSS
#'   \item Other original metadata
#' }
#' @importFrom GenomicRanges promoters distanceToNearest makeGRangesFromDataFrame mcols
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @export
#' @examples
#' # Annotate REs using nearest TSS and classify as Promoter or Distal
#' data(atac)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' retsi <- compute_spicey_index(atac=atac)
#' retsi_gene_nearest <- annotate_with_nearest(
#'   retsi = retsi,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db::org.Hs.eg.db
#' )
annotate_with_nearest <- function(retsi,
                                  txdb,
                                  annot_dbi,
                                  protein_coding_only = TRUE,
                                  keep_mito = FALSE,
                                  verbose = TRUE,
                                  add_tss_annotation = FALSE) {

  if (verbose) {
    message("Annotating regulatory elements to nearest gene...")
  }

  # Convert to GRanges if needed
  if (inherits(retsi, "data.frame") && !inherits(retsi, "GRanges")) {
    retsi <- GenomicRanges::makeGRangesFromDataFrame(retsi,
                                                     keep.extra.columns = TRUE)
  }

  # Remove ALT contigs if needed
  alt_chroms <- grep("_alt|random|fix|Un", seqlevels(retsi), value = TRUE)
  retsi <- GenomeInfoDb::keepSeqlevels(retsi,
                                       setdiff(seqlevels(retsi), alt_chroms),
                                       pruning.mode = "coarse")

  # Extract promoter regions
  proms <- get_promoters(txdb = txdb,
                         annot_dbi = annot_dbi,
                         keep_mito = keep_mito,
                         protein_coding_only = protein_coding_only,
                         verbose = verbose)

  # Find nearest TSS and compute distance
  nearest_hits <- GenomicRanges::distanceToNearest(retsi, proms)

  # Initialize annotation columns
  retsi$distanceToTSS <- NA_integer_
  retsi$nearestGeneSymbol <- NA_character_

  # Fill distance and gene symbol info for matched peaks
  retsi$distanceToTSS[queryHits(nearest_hits)] <- GenomicRanges::mcols(nearest_hits)$distance
  retsi$nearestGeneSymbol[queryHits(nearest_hits)] <- proms$gene_id[subjectHits(nearest_hits)]

  # Annotate as Promoter vs. Distal
  retsi$annotation <- ifelse(
    !is.na(retsi$distanceToTSS) & abs(retsi$distanceToTSS) <= 2000,
    "Promoter",
    "Distal"
  )

  # Optionally add direct TSS overlap annotation
  if (add_tss_annotation) {
    retsi <- annotate_tss(retsi,
                          txdb = txdb,
                          annot_dbi = annot_dbi)
  }

  # Remove names and convert to data.frame before returning
  retsi_df <- retsi |>
    data.frame()

  return(retsi_df)
}
