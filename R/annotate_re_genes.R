#' Link ATAC regions to putative target genes
#' @param annot_method Optional annotation method: \code{"nearest"} or \code{"coaccessibility"}.
#' @param links \code{GInteractions} object with co-accessibility links (required if \code{annot_method = "coaccessibility"}).
#' @param link_spicey_measures Logical; link RETSI and GETSI scores (default FALSE).
#' @param coaccess_cutoff_override Numeric; cutoff for co-accessibility clustering (default 0.25).
#' @param filter_promoter_distal Logical; filter to promoter-distal links (default TRUE).
#' @param filter_protein_coding Logical; restrict to protein-coding genes (default TRUE).
#' @param keep_mito Logical; keep mitochondrial chromosomes (default FALSE).
#' @param txdb \code{TxDb} object for genome annotation (required if annotation requested).
#' @param annot_dbi \code{AnnotationDbi} object for gene ID mapping (required if annotation requested).
#' @param add_tss_annotation Logical; annotate regulatory elements overlapping TSS (default FALSE).
#' @export
link_atac_to_genes <- function(granges=NULL, 
                               gene_id=NULL,
                               annot_method = NULL,
                               links = NULL,
                               link_spicey_measures = FALSE,
                               coaccess_cutoff_override = 0.25,
                               filter_promoter_distal = TRUE,
                               filter_protein_coding = TRUE,
                               txdb = NULL,
                               keep_mito = FALSE,
                               annot_dbi = NULL,
                               add_tss_annotation = FALSE) {
  
  
  if (is.null(txdb)) stop("Annotation requires a 'txdb' object.")
  if (is.null(annot_dbi)) stop("Annotation requires an 'annot_dbi' object.")
  if (is.null(retsi)) stop("ATAC data (RETSI) is required for region-to-gene annotation when 'annot_method' is specified.")
  
  retsi_annotated <- switch(
    annot_method,
    nearest = annotate_with_nearest(
      retsi = retsi,
      txdb = txdb,
      keep_mito = keep_mito,
      annot_dbi = annot_dbi,
      protein_coding_only = filter_protein_coding,
      verbose = verbose,
      add_tss_annotation = add_tss_annotation
    ),
    coaccessibility = {
      if (is.null(links)) stop("'links' must be provided for coaccessibility annotation.")
      annotate_with_coaccessibility(
        links = links,
        retsi = retsi,
        txdb = txdb,
        keep_mito = keep_mito,
        annot_dbi = annot_dbi,
        protein_coding_only = filter_protein_coding,
        verbose = verbose,
        coaccess_cutoff_override = coaccess_cutoff_override,
        filter_promoter_distal = filter_promoter_distal,
        add_tss_annotation = add_tss_annotation
      )
    },
    stop("Invalid 'annot_method'. Choose 'nearest' or 'coaccessibility'.")
  )
  
  if (link_spicey_measures) {
    if (is.null(getsi)) stop("RNA data must be provided to link RETSI and GETSI.")
    combined <- link_spicey(
      retsi_annotated = retsi_annotated,
      getsi = getsi,
      method = annot_method)
    message("SPICEY pipeline successfully completed")
    return(combined)
  }
}


#' Annotate regulatory elements overlapping transcription start sites (TSS)
#'
#' Identifies regulatory elements that overlap precisely defined transcription
#' start sites (TSS) and assigns the corresponding gene symbols to these elements.
#' @param txdb A `TxDb` object (from \pkg{GenomicFeatures}), a `GRanges` object,
#' or any object that supports `GenomeInfoDb::seqlevels()`.
#' For example, you can create a TxDb using: `GenomicFeatures::makeTxDbFromGFF()`
#' or use a prebuilt TxDb such as `TxDb.Hsapiens.UCSC.hg38.knownGene`.
#' @param re A \code{GRanges} or \code{data.frame} with regulatory regions to be annotated.
#'           If a \code{data.frame}, it must be convertible to \code{GRanges} with genomic coordinates.
#' @param annot_dbi An `AnnotationDbi` object (e.g., \pkg{org.Hs.eg.db})
#' for mapping gene IDs to gene symbols and types.
#' @return A \code{data.frame} containing the input regulatory elements with added columns:
#' \describe{
#'   \item{\code{in_TSS}}{Logical indicating whether the element overlaps a TSS.}
#'   \item{\code{TSS_gene}}{Gene symbol of the overlapping TSS, or \code{NA} if none.}}
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

#' Annotate co-accessible links with CCAN membership
#'
#' Processes co-accessibility links to assign CCAN (co-accessibility network)
#' membership to each peak. Optionally filters links to retain only promoter-distal
#' interactions, and converts the result to a \code{GInteractions} object annotated
#' with co-accessibility scores and CCAN metadata.
#' @param links A \code{data.frame} of co-accessibility links with at least three columns:
#' \describe{
#'   \item{\code{Peak1}}{Character vector specifying the first peak in the
#'   co-accessibility link. Each peak is typically represented as a
#'   genomic coordinate string, e.g., "chr1-123456-123789".}
#'   \item{\code{Peak2}}{Character vector specifying the second peak in the
#'   co-accessibility link, in the same format as \code{Peak1}.}
#'   \item{\code{coaccess}}{Numeric vector representing the co-accessibility
#'   score between \code{Peak1} and \code{Peak2}.
#'   This score quantifies the strength of coordinated accessibility between the
#'   two peaks, where higher values indicate stronger putative regulatory interaction.}}
#' @param coaccess_cutoff_override Numeric; co-accessibility score cutoff used
#'   for defining CCAN membership clusters. Peaks connected by edges with co-accessibility
#'   scores above this threshold will be clustered together (default: \code{0.25}).
#' @param tolerance_digits Integer; number of decimal digits to round co-accessibility
#'   scores when detecting CCAN clusters. This controls the precision of clustering (default: \code{2}).
#' @param filter_promoter_distal Logical; whether to filter links to keep only those connecting
#'   promoter regions to distal elements (default: \code{TRUE}). Promoter regions are defined
#'   based on the provided \code{txdb} annotation.
#' @param txdb A \code{TxDb} object (from \pkg{GenomicFeatures}) representing genome annotation.
#'   This is used to identify promoter regions when filtering promoter-distal links. Can be created using
#'   \code{GenomicFeatures::makeTxDbFromGFF()} or a prebuilt TxDb package like \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' @return A \code{GInteractions} object where each row corresponds to a co-accessible peak pair,
#'   with the following metadata columns:
#'   \describe{
#'     \item{\code{coaccess}}{Co-accessibility score between \code{Peak1} and \code{Peak2}.}
#'     \item{\code{CCAN1}}{Identifier of the CCAN (co-accessibility network cluster)
#'     that \code{Peak1} belongs to. CCANs group highly connected peaks based on co-accessibility.}
#'     \item{\code{CCAN2}}{Identifier of the CCAN cluster that \code{Peak2} belongs to.}}
#' @importFrom dplyr left_join filter rename
#' @importFrom regioneR toGRanges
#' @importFrom ChIPseeker annotatePeak
#' @importFrom InteractionSet GInteractions
#' @importFrom Signac StringToGRanges
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




#' Extract CCAN membership for genomic regions from co-accessibile Links
#'
#' Given a set of co-accessibility links with CCAN (co-accessibility network) annotations,
#' this function finds which CCANs each genomic region belongs to, or conversely, which
#' regions belong to each CCAN. This helps to identify the grouping of regulatory elements
#' within shared co-accessibility networks.
#' @param links A \code{GInteractions} object containing co-accessibility links.
#'              Must have metadata columns \code{CCAN1} and \code{CCAN2}
#'              representing CCAN IDs for each anchor in the interaction.
#' @param gr A \code{GRanges} object representing genomic regions to query.
#' @param name_column A \code{character} string specifying the metadata column name in
#'                    \code{gr} to use as the region identifier. Default is \code{"region"}.
#' @param split A \code{character} indicating the desired output split:
#'              \itemize{
#'                \item \code{"ccan"}: returns a named list with CCAN IDs as names,
#'                       each containing a vector of region identifiers.
#'                \item \code{"name"}: returns a named list with region identifiers as names,
#'                       each containing a vector of CCAN IDs.
#'              }
#'              Default is \code{c("ccan", "name")}, with \code{"ccan"} selected by default.
#' @return A named \code{list} mapping CCANs to region names or vice versa, depending on
#'         the \code{split} argument.
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
#'
#' This function annotates regulatory elements by linking them to gene promoters
#' through shared CCAN (co-accessibility network) membership derived from co-accessibile links.
#' Each regulatory region is assigned the gene IDs of promoters that belong to the same CCAN,
#' thus identifying putative gene targets based on chromatin co-accessibility.
#' @param links A \code{GInteractions} object containing co-accessibility links.
#'              Must have metadata columns \code{CCAN1} and \code{CCAN2}
#'              representing CCAN IDs for each anchor in the interaction.
#' @param re A \code{GRanges} object of regulatory elements with a metadata column
#'           named \code{"region"} that identifies each region.
#' @param proms A \code{GRanges} object of gene promoters
#' @return A \code{GRanges} object identical to \code{re} but with an added metadata
#'         column \code{gene_coacc} that contains a list of linked gene IDs per RE.
get_targets_links <- function(links, re, proms) {
  res <- get_ccan(links, re, name_column = "region", split = "name")
  re$CCAN <- res[re$region]
  genes <- get_ccan(links, proms, name_column = "gene_id", split = "ccan")
  GenomicRanges::mcols(re)$gene_coacc <- lapply(re$CCAN, function(x) {
    if (!is.null(x)) unlist(genes[as.character(x)]) else NA
  })
  return(re)
}




#' Annotate regulatory regions to their gene targets using co-accessible links
#'
#' Annotate a set of regulatory elements with their putative target genes
#' using co-accessibility links annotated with CCAN (co-accessibility network) metadata.
#' The function extracts promoter regions from a given \code{TxDb} annotation,
#' optionally filters for protein-coding genes, and annotates REs based on
#' shared CCAN membership with promoters. Additionally, it computes the distance
#' of each RE to the nearest TSS and can add direct TSS overlap annotations.
#' @param links A \code{GInteractions} object containing co-accessibility links.
#'   Must have metadata columns \code{CCAN1} and \code{CCAN2}
#'   representing CCAN IDs for each anchor in the interaction.
#' @param retsi A \code{GRanges} or \code{data.frame} of regulatory elements to annotate.
#'   If a \code{data.frame}, it must have columns sufficient to construct
#'   a \code{GRanges} object, including \code{seqnames}, \code{start}, and \code{end}.
#' @param txdb A \code{TxDb} object (from \pkg{GenomicFeatures}) representing genome annotation.
#'   Can be created using \code{GenomicFeatures::makeTxDbFromGFF()} or a prebuilt
#'   TxDb package such as \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' @param annot_dbi An \code{AnnotationDbi} object (e.g., \pkg{org.Hs.eg.db})
#'   for mapping gene IDs to gene symbols and types.
#' @param protein_coding_only Logical, default \code{TRUE}.
#'   If \code{TRUE}, restricts to protein-coding genes based on the \code{GENETYPE} annotation in \code{annot_dbi}.
#' @param keep_mito Logical, default \code{FALSE}. Whether to keep mitochondrial chromosomes.
#' @param verbose Logical, default \code{TRUE}. If \code{TRUE}, prints informative messages.
#' @param coaccess_cutoff_override Numeric; co-accessibility score cutoff used
#'   for defining CCAN membership clusters. Default is \code{0.25}.
#' @param filter_promoter_distal Logical, default \code{TRUE}.
#'   Whether to filter links to keep only those connecting promoter regions to distal elements.
#' @param add_tss_annotation Logical, default \code{TRUE}.
#'   If \code{TRUE}, annotate REs that overlap TSS regions directly.
#' @return A \code{data.frame} of regulatory elements annotated with:
#' \itemize{
#'   \item Distance to nearest TSS and classification as "Promoter" or "Distal".
#'   \item Putative target genes linked by co-accessibility.
#'   \item Optional direct TSS overlap and associated gene symbols.}
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
#'   coaccess_cutoff_override = 0.25)
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
  if (inherits(retsi, "data.frame") && !inherits(retsi, "GRanges")) {
    retsi <- GenomicRanges::makeGRangesFromDataFrame(retsi,
                                                     keep.extra.columns = TRUE)
  }
  alt_chroms <- grep("_alt|random|fix|Un", seqlevels(retsi), value = TRUE)
  retsi <- GenomeInfoDb::keepSeqlevels(retsi,
                                       setdiff(seqlevels(retsi), alt_chroms),
                                       pruning.mode = "coarse")
proms <- get_promoters(txdb = txdb,
                         annot_dbi = annot_dbi,
                         keep_mito = keep_mito,
                         protein_coding_only = protein_coding_only,
                         verbose = verbose)
  nearest_hits <- GenomicRanges::distanceToNearest(retsi, proms)
  retsi$distanceToTSS <- NA_integer_
  retsi$distanceToTSS[queryHits(nearest_hits)] <- GenomicRanges::mcols(nearest_hits)$distance
  retsi$annotation <- ifelse(
    !is.na(retsi$distanceToTSS) &
    abs(retsi$distanceToTSS) <= 2000, "Promoter","Distal")

  annotated_links <- annotate_links_with_ccans(
    links = links,
    coaccess_cutoff_override = coaccess_cutoff_override,
    filter_promoter_distal = filter_promoter_distal,
    txdb = txdb)

  re_annotated <- get_targets_links(
    links = annotated_links,
    re = retsi,
    proms = proms)

  re_unnested <- tidyr::unnest(as.data.frame(re_annotated),
                               cols = "gene_coacc") |>
    dplyr::select(-c(CCAN)) |>
    data.frame()

  if (add_tss_annotation) {
    re_unnested <- annotate_tss(re_unnested,
                                txdb = txdb,
                                annot_dbi = annot_dbi)
  }
  return(re_unnested)
}






#' Annotate regulatory regions to their gene targets using distance to the nearest transcription start site (TSS)
#'
#' Annotate a set of regulatory elements with their putative target genes
#' based on distance to the nearest TSS, classifying as "Promoter" (<2kb) or "Distal",
#' and adding the nearest gene symbol using Biomart-annotated genes.
#' Optionally, annotate direct overlap with TSS using 1-bp TSS coordinates.
#' @param retsi A \code{GRanges} or \code{data.frame} of regulatory elements to annotate.
#'   If a \code{data.frame}, it must have columns sufficient to construct
#'   a \code{GRanges} object, including \code{seqnames}, \code{start}, and \code{end}.
#' @param txdb A \code{TxDb} object (from \pkg{GenomicFeatures}) representing genome annotation.
#'   Can be created using \code{GenomicFeatures::makeTxDbFromGFF()} or a prebuilt
#'   TxDb package such as \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' @param annot_dbi An \code{AnnotationDbi} object (e.g., \pkg{org.Hs.eg.db})
#'   for mapping gene IDs to gene symbols and types.
#' @param protein_coding_only Logical, default \code{TRUE}.
#'   If \code{TRUE}, restricts to protein-coding genes based on the \code{GENETYPE} annotation in \code{annot_dbi}.
#' @param keep_mito Logical, default \code{FALSE}. Whether to keep mitochondrial chromosomes.
#' @param verbose Logical, default \code{TRUE}. If \code{TRUE}, prints informative messages.
#' @param add_tss_annotation Logical, default \code{TRUE}.
#'   If \code{TRUE}, annotate REs overlapping TSS directly.
#' @return A \code{data.frame} containing:
#' \itemize{
#'   \item \code{distanceToTSS}: Distance to nearest TSS.
#'   \item \code{annotation}: "Promoter" (< 2kb) or "Distal".
#'   \item \code{nearestGeneSymbol}: Symbol of the nearest gene.
#'   \item Optional \code{in_TSS} and \code{TSS_gene} if \code{add_tss_annotation = TRUE}.
#'   \item Other original metadata columns.}
#' @importFrom GenomicRanges promoters distanceToNearest makeGRangesFromDataFrame mcols
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @export
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#' data(atac)
#' retsi <- compute_spicey_index(atac = atac)
#' retsi_gene_nearest <- annotate_with_nearest(
#'   retsi = retsi,
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   annot_dbi = org.Hs.eg.db)
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
  if (inherits(retsi, "data.frame") && !inherits(retsi, "GRanges")) {
    retsi <- GenomicRanges::makeGRangesFromDataFrame(retsi,
                                                     keep.extra.columns = TRUE)
  }
  alt_chroms <- grep("_alt|random|fix|Un", seqlevels(retsi), value = TRUE)
  retsi <- GenomeInfoDb::keepSeqlevels(retsi,
                                       setdiff(seqlevels(retsi), alt_chroms),
                                       pruning.mode = "coarse")

  proms <- get_promoters(txdb = txdb,
                         annot_dbi = annot_dbi,
                         keep_mito = keep_mito,
                         protein_coding_only = protein_coding_only,
                         verbose = verbose)

  nearest_hits <- GenomicRanges::distanceToNearest(retsi, proms)
  retsi$distanceToTSS <- NA_integer_
  retsi$nearestGeneSymbol <- NA_character_
  retsi$distanceToTSS[queryHits(nearest_hits)] <- GenomicRanges::mcols(nearest_hits)$distance
  retsi$nearestGeneSymbol[queryHits(nearest_hits)] <- proms$gene_id[subjectHits(nearest_hits)]

  retsi$annotation <- ifelse(
    !is.na(retsi$distanceToTSS) &
    abs(retsi$distanceToTSS) <= 2000, "Promoter","Distal")

  if (add_tss_annotation) {
    retsi <- annotate_tss(retsi,
                          txdb = txdb,
                          annot_dbi = annot_dbi)
  }
  retsi_df <- retsi |>
    data.frame()
  return(retsi_df)
}
