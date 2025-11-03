#' Helper function to generate mutated barcodes for testing and examples
#'
#' @param barcodes A \code{\link{DNAStringSet}} containing
#' the barcodes to be mutated
#' @param mismatches Integer vector, number of errors to make
#' @param times Integer vector the same length as \code{mismatches}: The number of
#' sequences having the corresponding number of errors
#'
#' @returns A list with the following elements
#' \itemize{
#' \item \code{barcodes}: A \code{\link{DNAStringSet}} with the mutated barcodes.
#' All barcodes are assumed to have the same width
#' \item \code{mismatches}: An integer vector with
#' the same length as \code{barcodes} showing the number
#' of errors introduced in each barcode. Remaining barcodes will
#' be returned without errors
#' }
#' @importFrom Biostrings width
#' @importFrom utils combn
#'
#' @noRd
mutate_barcodes <- function(barcodes, mismatches, times) {
    widths <- width(barcodes)
    unique_width <- unique(widths)
    assert_that(length(unique_width) <= 1L, msg = glue("All barcodes must have the same width"))
    assert_that(length(mismatches) == length(times), msg = "Arguments mismatches and times must have the same length")
    n_barcodes <- length(barcodes)
    assert_that(
        sum(times) <= n_barcodes,
        msg = "Total number of mutated barcodes must be less than \
        or equal to the number of provided barcodes"
    )
    barcode_mismatches <- rep(0L, n_barcodes)
    # A boolean flag vector showing whether an index has been sampled from
    idx_available <- rep(TRUE, n_barcodes)
    sampling_bases <- rep_len(Biostrings::DNA_BASES, Biostrings::DNA_BASES %>%
        length() %>%
        add(. - 1L))
    barcode_matrix <- barcodes %>%
        as.character() %>%
        strsplit(split = "") %>%
        unlist() %>%
        matrix(nrow = unique_width, ncol = n_barcodes)
    index_matrix <- match(barcode_matrix, Biostrings::DNA_BASES) %>%
        matrix(ncol = n_barcodes, dimnames = dimnames(barcode_matrix))
    # I think this is a valid use case for explicit loops in R
    for (i in seq_along(mismatches)) {
        this_mismatches <- mismatches[i]
        idxs_with_this_mismatches <- sample(which(idx_available), size = times[i])
        idx_available[idxs_with_this_mismatches] <- FALSE
        barcode_mismatches[idxs_with_this_mismatches] <- this_mismatches
        possible_error_positions <- combn(seq_len(unique_width), this_mismatches)
        selected_error_positions <- possible_error_positions[, x = sample(possible_error_positions %>%
            ncol(), size = times[i], replace = TRUE)]
        # We add a random offset of length 3 to the indicies of the DNA bases
        base_offsets <- sample(x = Biostrings::DNA_BASES %>%
            length() %>%
            subtract(1L) %>%
            seq_len(), size = this_mismatches * times[i], replace = TRUE)
        subscript_matrix <- cbind(selected_error_positions %>%
            as.vector(), rep(idxs_with_this_mismatches, each = this_mismatches))
        index_matrix[subscript_matrix] <- index_matrix[subscript_matrix] + base_offsets
    }
    # We could paste everything together, but we choose to optimize the process such
    # that only reads which are mutated are decoded and pasted
    mutated_barcodes <- which(!idx_available)
    barcodes[mutated_barcodes] <- index_matrix[, mutated_barcodes] %>%
        {
            array(sampling_bases[.], dim = dim(.))
        } %>%
        apply(MARGIN = 2L, FUN = paste0, collapse = "")
    names(barcode_mismatches) <- names(barcodes)
    list(barcodes = barcodes, mismatches = barcode_mismatches)
}
