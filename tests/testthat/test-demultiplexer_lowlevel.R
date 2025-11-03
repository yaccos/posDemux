library(Biostrings)

barcodes <- DNAStringSet(c(a = "ATCG", b = "GGGA", c = "CCTA"))
segments <- DNAStringSet(c(w = "TAAC", x = "ATCG", y = "GCGA", z = "CCAT"))
assigned_barcodes <- c(w = "a", x = "a", y = "b", z = "c")
mismatches <- c(w = 4L, x = 0L, y = 1L, z = 2L)
test_that("Hamming matching works", {
    demux_res <- posDemux:::hamming_match(
        segments, names(segments), barcodes, names(barcodes),
        4L
    )
    expect_equal(demux_res$assigned_barcodes, assigned_barcodes)
    expect_equal(demux_res$mismatches, mismatches)
})

test_that("Hamming matching gives expected output on unnamed segments", {
    segments <- unname(segments)
    demux_res <- posDemux:::hamming_match(
        segments, names(segments), barcodes, names(barcodes),
        4L
    )
    expect_equal(demux_res$assigned_barcodes, unname(assigned_barcodes))
    expect_equal(demux_res$mismatches, unname(mismatches))
})
