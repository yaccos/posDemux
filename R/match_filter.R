#' Filter demultiplexed reads
#'
#' @description
#' Filters the demultiplexed reads from \code{\link{combinatorial_demultiplex}}
#' such that any read exceeding the number of allowed mismatches for any
#' of the barcodes is removed. The function gives diagnostic information
#' on the number of reads removed per barcode and the total number
#' of reads removed.
#'
#' @param demultiplex_res Unprocessed output from
#'  \code{\link{combinatorial_demultiplex}}
#' @param allowed_mismatches Integer vector of length one or the same length
#' as the number of barcode segments, the maximum Hamming distance from including
#' a read in the output.
#' @returns A list with the following elements:
#' \itemize{
#' \item \code{demultiplex_res}: The contents of \code{demultiplex_res} with
#' the sequences filtered.
#' \item \code{retained}: Logical vector with the same length as
#'  the number of reads in the input. \code{TRUE} if the corresponding read
#'  is retained.
#'  Useful for future filtering of paired-end reads.
#' \item \code{summary_res}: Result of \code{\link{create_summary_res}} called
#' on the results of filtering.
#' }
#' @details
#' The value of \code{n_removed} does not in general equal the sum of
#' \code{n_removed_per_barcode} since a read can have too many mismatches
#' with multiple barcodes.
#'
#' @importFrom magrittr set_names %<>%
#' @importFrom purrr imap
#' @export
#'
filter_demultiplex_res <- function(demultiplex_res, allowed_mismatches) {
  barcodes <- demultiplex_res$barcodes
  n_barcode_sets <- length(barcodes)
  assert_that(
    length(allowed_mismatches) == 1L ||
      n_barcode_sets == length(allowed_mismatches),
    msg = "allowed_mismatches does not match the number of barcode sets"
  )
  if (length(allowed_mismatches) == 1L && n_barcode_sets > 1L) {
    allowed_mismatches <- rep(allowed_mismatches, n_barcode_sets)
  }
  names(allowed_mismatches) <- names(demultiplex_res$barcodes)
  mismatches <- demultiplex_res$mismatches
  raw_filter_res <- filter_sequences(demultiplex_res, allowed_mismatches)
  retained_sequences <- raw_filter_res$retained_sequences
  res <- raw_filter_res$res
  summary_res <- create_summary_res(retained_sequences,
                                    barcodes,
                                    res$assigned_barcodes,
                                    allowed_mismatches,
                                    mismatches)
  
  list(
    demultiplex_res = res,
    retained = retained_sequences,
    summary_res = summary_res
  )
}

filter_sequences <- function(demultiplex_res, allowed_mismatches) {
  mismatches <- demultiplex_res$mismatches
  mismatches_above_threshold <- sweep(mismatches, 2L, allowed_mismatches, FUN = `>`)
  retained_sequences <- rowSums(mismatches_above_threshold) == 0L
  demultiplex_res$assigned_barcodes %<>% extract(retained_sequences, )
  demultiplex_res$mismatches %<>% extract(retained_sequences, )
  demultiplex_res$payload %<>% map(. %>% extract(retained_sequences))
  list(res = demultiplex_res, retained_sequences = retained_sequences)
}
