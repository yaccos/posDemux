library(Biostrings)
library(purrr)
reads <- DNAStringSet()
barcode_files <- system.file("extdata/PETRI-seq_barcodes", c(
    bc1 = "bc1.fa", bc2 = "bc2.fa",
    bc3 = "bc3.fa"
), package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
output_table_file <- tempfile(fileext = ".txt")
barcode_index <- map(barcode_files, readDNAStringSet)
barcodes <- barcode_index[c("bc3", "bc2", "bc1")]
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
demultiplex_res <- combinatorial_demultiplex(reads, barcodes,
                                             sequence_annotation, segment_lengths)
expected_dimnames <- list(NULL, c(c("bc3", "bc2", "bc1")))
expected_assigned_barcodes <- matrix(
    character(), nrow = 0L, ncol = 3L, dimnames = expected_dimnames)
expected_mismatches <- matrix(
    integer(), nrow = 0L, ncol = 3L, dimnames = expected_dimnames)
expected_payload <- list(UMI = DNAStringSet())
expected_demultiplex_res <- list(
    assigned_barcodes = expected_assigned_barcodes,
    mismatches = expected_mismatches,
    payload = expected_payload,
    barcodes = barcodes
    )

test_that("Empty stringsets are demultiplexed as expected", {
    expect_equal(demultiplex_res, expected_demultiplex_res)
})
