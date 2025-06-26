#' Annotate and Filter Gene Expression Differential Data, Return GRanges
#'
#' This function takes a named list of gene expression differential analysis results
#' (typically one data frame per cell type). It maps gene symbols to Ensembl IDs,
#' joins with genomic coordinates from Biomart, filters to protein-coding genes
#' on canonical chromosomes, converts to a GRanges object, and adds a `cell_type`
#' metadata column based on the list name.
#'
#' @param rna_da A named list of data frames, where each element corresponds to a cell type
#'        and contains differential gene expression results. Row names should be gene symbols.
#' @return A single GRanges object combining all annotated genes, with metadata columns
#'         including differential statistics and a `cell_type` column.
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr left_join filter select mutate everything
#' @importFrom GenomicRanges makeGRangesFromDataFrame sort mcols
#' @importFrom purrr imap reduce
#'
#' @examples
#' \dontrun{
#' combined_gr <- combine_gex_da(rna_da)
#' }
#'
#' @export
combine_gex_da <- function(rna_da) {
  # Annotate each element of the list
  gr_list_annot <- lapply(names(rna_da), function(cell_type) {
    df <- rna_da[[cell_type]]

    # Ensure gene symbols are a column
    df$symbol <- rownames(df)

    # Add cell_type column if missing
    if (!"cell_type" %in% colnames(df)) {
      df$cell_type <- cell_type
    }

    # Map SYMBOL -> ENSEMBL using org.Hs.eg.db
    df$ensembl_id <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = df$symbol,
      column = "ENSEMBL",
      keytype = "SYMBOL",
      multiVals = "first"
    )

    # Join with Biomart gene info
    df_annot <- dplyr::left_join(
      df,
      biomart_genes()$df |> dplyr::select(
        chromosome_name,
        start_position,
        end_position,
        strand,
        ensembl_gene_id,
        gene_biotype
      ),
      by = c("ensembl_id" = "ensembl_gene_id")
    )

    if (nrow(df_annot) == 0) {
      message("No annotation found for cell type: ", cell_type)
      return(NULL)
    }

    # Filter and convert to GRanges
    df_annot |>
      dplyr::filter(
        !is.na(chromosome_name),
        !is.na(start_position),
        !is.na(end_position),
        gene_biotype == "protein_coding",
        chromosome_name %in% c(as.character(1:22), "X", "Y")
      ) |>
      dplyr::mutate(chromosome_name = paste0("chr", chromosome_name)) |>
      GenomicRanges::makeGRangesFromDataFrame(
        seqnames.field = "chromosome_name",
        start.field = "start_position",
        end.field = "end_position",
        strand.field = "strand",
        keep.extra.columns = TRUE
      ) |>
      GenomicRanges::sort()
  })

  # Name each GRanges with cell type
  names(gr_list_annot) <- names(rna_da)

  # Drop NULLs (missing annotations)
  gr_list_annot <- Filter(Negate(is.null), gr_list_annot)

  # Add `cell_type` metadata and merge into one GRanges object
  combined_gr <- purrr::imap(gr_list_annot, ~ {
    GenomicRanges::mcols(.x)$cell_type <- .y
    .x
  }) |>
    purrr::reduce(c)

  return(combined_gr)
}





#' Compute GETSI scores on RNA data
#'
#' @param gr GRanges object with RNA differential expression data (avg_log2FC, p_val, cell_type, symbol)
#' @return GRanges with computed GETSI scores (GETSI column)
#' @importFrom regioneR toGRanges
#' @import dplyr
#' @export
getsi <- function(rna_da) {
  rna_da_gr <- combine_gex_da(rna_da)
  getsi_gr <- data.frame(rna_da_gr) |>
    mutate(
      region = paste0(seqnames, ":", start, "-", end),
      avg_FC = 2^(avg_log2FC),
      p_val_adj = p.adjust(p_val, method = "fdr")
    ) |>
    group_by(symbol) |>
    mutate(
      max_log2FC = max(avg_FC, na.rm = TRUE)
      # max_log2FC = ifelse(max_log2FC == 0, 1e-6, max_log2FC)
    ) |>
    ungroup() |>
    group_by(cell_type) |>
    mutate(
      p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm=TRUE), p_val_adj),
      weight = scales::rescale(-log10(p_val_adj), to=c(0, 1))
    ) |>
    group_by(cell_type, symbol) |>
    mutate(
      norm_log2FC = avg_FC / max_log2FC,
      GETSI = (norm_log2FC * weight)
    ) |>
    ungroup() |>
    dplyr::select(-c(width, strand)) |>
    data.frame() |>
    regioneR::toGRanges()
  return(getsi_gr)
}



#' Compute normalized entropy of GETSI scores
#'
#' @param gr GRanges with GETSI scores and symbol column
#' @return Data frame with symbol and norm_entropy columns
#' @import dplyr
#' @export
entropy_getsi <- function(getsi_gr) {
  df <- data.frame(getsi_gr)
  n_celltypes <- length(unique(df$cell_type))

  df <- df |>
    group_by(symbol) |>
    mutate(
      GETSI_prob = GETSI / sum(GETSI, na.rm = TRUE),
      entropy_component = -GETSI_prob * log2(GETSI_prob)
    ) |>
    summarise(
      entropy = sum(entropy_component, na.rm = TRUE),
      # norm_entropy = entropy / log2(n_celltypes),
      norm_entropy = 1 - exp(-entropy),
      .groups = "drop"
    )
  return(df)
}


#' Compute GETSI scores and entropy from RNA-seq data
#'
#' Computes gene expression tissue specificity (GETSI) scores and
#' entropy from RNA-seq data. The input should be a named list of
#' GRanges objects, each containing differential expression data for
#' one cell type.
#'
#' @param rna_da A named list of GRanges objects, typically one per cell type,
#'   containing RNA differential expression data with columns such as
#'   avg_log2FC, p_val, symbol, etc.
#'
#' @return A GRanges object with GETSI scores and entropy values in the
#'   metadata columns.
#'
#' @export
spicey_getsi <- function(rna_da) {
  message("→ Combining, annotating, and computing GETSI scores...")
  getsi_gr <- getsi(rna_da)

  message("→ Computing entropy on GETSI scores...")
  entropy_df <- entropy_getsi(getsi_gr)

  # Join entropy info with getsi results
  getsi_final <- as.data.frame(getsi_gr) |>
    dplyr::left_join(entropy_df, by = "symbol") |>
    dplyr::select(-c(width, strand, avg_FC,
                     max_log2FC, norm_log2FC,
                     weight, entropy)) |>
    regioneR::toGRanges()

  # Return combined GRanges with GETSI and entropy columns
  return(getsi_final)
}
