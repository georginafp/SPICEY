# ...............
source("R/utils.R")

# ATAC ----
merged_df <- readRDS("/homes/users/gfuentes/scratch/projects/SPICEY/data/ATAC_DAR.rds")
merged_df <- unlist(GRangesList(merged_df))

retsi_gr <- merged_df %>%
  data.frame() %>%
  dplyr::select(c(seqnames, start, end, everything())) %>%
  dplyr::select(-c(width, strand)) %>%
  regioneR::toGRanges()

anno <- ChIPseeker::annotatePeak(retsi_gr,
                                 TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                                 verbose=FALSE)

retsi_gr$distanceToTSS <- data.frame(anno)$distanceToTSS
retsi_gr$annotation <- "Distal"
retsi_gr$annotation[abs(retsi_gr$distanceToTSS)<=2e3] <- "Promoter"

genes <- biomart_genes()$gr

retsi_gr <- retsi_gr %>%
  plyranges::mutate(nearestGeneSymbol = genes$external_gene_name[nearest(retsi_gr, promoters(genes, 1, 0))])
names(retsi_gr) <- NULL
saveRDS(retsi_gr, "data/FINAL_ATAC.rds")


# RNA ----
gr_list <- readRDS("/homes/users/gfuentes/scratch/projects/SPICEY/data/RNA_DEG.rds")


gr_df_list <- map(gr_list, ~ as.data.frame(.x))
gr_list <- lapply(names(gr_df_list), function(name) {

  df <- gr_df_list[[name]]
  df$symbol <- rownames(df)
  df$ensembl_id <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                         keys=df$symbol,
                                         column="ENSEMBL",
                                         keytype="SYMBOL",
                                         multiVals="first")

  df_annot <- dplyr::left_join(df,
                               biomart_genes()$df %>% dplyr::select(chromosome_name,
                                                                    start_position,,
                                                                    end_position,
                                                                    strand,
                                                                    ensembl_gene_id,
                                                                    gene_biotype),
                               by=c(ensembl_id="ensembl_gene_id"))

  if (nrow(df_annot) == 0) {
    message("Missing entry annotation for ", name)
    return(NULL)
  }

  df_annot %>%
    dplyr::filter(!is.na(chromosome_name),
                  !is.na(start_position),
                  !is.na(end_position),
                  gene_biotype == "protein_coding") %>%
    dplyr::filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
    dplyr::select(chromosome_name,
                  start_position,
                  end_position,
                  strand,
                  everything()) %>%
    dplyr::mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
    makeGRangesFromDataFrame(seqnames.field=c("seqnames", "seqname",
                                              "chromosome", "chrom",
                                              "chr", "chromosome_name",
                                              "seqid"),
                             start.field=c("start", "start_position"),
                             end.field=c("end", "stop", "end_position"),
                             strand.field="strand",
                             keep.extra.columns = T) %>%
    sort()
})


# Rename list with cell types
names(gr_list) <- names(gr_df_list)
gr_list <- Filter(Negate(is.null), gr_list) # Remove NULLs


# Bind together in one GRanges
retsi_gr <- imap(gr_list, ~ {
  mcols(.x)$cell_type <- .y  # .x = GRanges, .y = name
  .x
}) %>%
  purrr::reduce(c) #214896

saveRDS(retsi_gr, "data/FINAL_RNA.rds")

