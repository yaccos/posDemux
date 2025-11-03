library(posDemux)
library(magrittr)
library(Biostrings)
set.seed(5203)

width <- 10L
n_barcodes <- length(LETTERS)
barcodes_matrix <- sample(DNA_BASES, size = width * n_barcodes, replace = TRUE) %>%
    matrix(ncol = n_barcodes) %>%
    set_colnames(LETTERS)

barcodes_stringset <- barcodes_matrix %>%
    apply(2L, FUN = paste0, collapse = "") %>%
    DNAStringSet()

# Up to three mismatches
mismatches <- seq_len(3L)
times <- c(10L, 6L, 3L)


mutated_barcodes <- posDemux:::mutate_barcodes(
    barcodes = barcodes_stringset, mismatches = mismatches,
    times = times
)
stringset_to_matrix <- function(stringset) {
    stringset %>%
        as.character() %>%
        strsplit(split = "") %>%
        unlist() %>%
        matrix(nrow = width, ncol = n_barcodes) %>%
        set_colnames(names(stringset))
}

hamming_distance <- function(x, y) {
    ((x %>%
        stringset_to_matrix()) != (y %>%
        stringset_to_matrix())) %>%
        colSums()
}


hamming_mismatches <- hamming_distance(barcodes_stringset, mutated_barcodes$barcodes)


test_that("Barcode mutation provides the reported number of mutations", {
    expect_equal(hamming_mismatches, mutated_barcodes$mismatches)
})

test_that("Barcode mutation gives the correct distribution of mutations", {
    mismatch_sum <- sum(hamming_mismatches)
    expect_equal(mismatch_sum, (mismatches * times) %>%
        sum())
    n_unmutated_barcodes <- length(barcodes_stringset) - sum(times)
    expect_equal((hamming_mismatches == 0) %>%
        sum(), n_unmutated_barcodes)
    expect_equal(outer(mismatches, hamming_mismatches, equals) %>%
        rowSums(), times)
})
