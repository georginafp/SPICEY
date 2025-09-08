data(rna)
data(atac)
data(cicero_links)

test_that("SPICEY integrates RNA, ATAC, and annotation", {
  peaks <- unique(unlist(atac)[, c("region_id")])
  annotation_coacc <- annotate_with_coaccessibility(
      peaks = peaks,
      txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
      links_df = cicero_links,
      annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
      protein_coding_only = TRUE,
      verbose = TRUE,
      add_tss_annotation = FALSE,
      upstream = 2000,
      downstream = 2000)

  res <- SPICEY(rna = rna, atac = atac, annotation = annotation_coacc)

  expect_type(res, "list")
  expect_true(all(c("RETSI", "GETSI", "linked") %in% names(res)))
})
