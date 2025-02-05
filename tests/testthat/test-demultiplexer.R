library(Biostrings)
test_that("Hamming matching works", {
  barcodes <- DNAStringSet(c(a="ATCG", b="GGGA", c="CCTA"))
  segments <- DNAStringSet(c(a="ATCG", b="GCGA", c="CCAT"))
  demux_res <- posDemux:::hamming_match(segments, names(segments),
                                        barcodes, names(barcodes), 4L)
  expect_equal(demux_res$assigned_barcode, c(a="a", b="b", c="c"))
  expect_equal(demux_res$mismatches, c(a=0L, b=1L, c=2L))
})
