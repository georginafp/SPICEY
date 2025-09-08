data(atac)

test_that("annotate_with_nearest works with example peaks", {
  peaks <- unique(unlist(atac)[, "region_id"])
  ann <- annotate_with_nearest(
    peaks = peaks,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    annot_dbi = org.Hs.eg.db::org.Hs.eg.db,
    protein_coding_only = TRUE,
    verbose = TRUE,
    add_tss_annotation = FALSE,
    upstream = 2000,
    downstream = 2000)

  expect_s3_class(ann, "data.frame")
  expect_true(all(c("region_id", "gene_id", "distanceToTSS", "annotation") %in% colnames(ann)))
})
