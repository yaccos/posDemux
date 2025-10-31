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
#'  \code{\link{combinatorial_demultiplex}}.
#' @param allowed_mismatches Integer vector of length one or the same length
#' as the number of barcode segments; the threshold Hamming distance. All reads
#' having a number of mismatches above this number in any of the barcodes will
#' be filtered away.
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
#' @importFrom purrr imap
#' @example inst/examples/match_filter-examples.R
#' @seealso [create_summary_res()]
#' @export
#'
filter_demultiplex_res <- function(demultiplex_res, allowed_mismatches) {
  barcodes <- demultiplex_res$barcodes
  allowed_mismatches <- validate_allowed_mismatches(allowed_mismatches, barcodes)
  mismatches <- demultiplex_res$mismatches
  raw_filter_res <- filter_sequences(demultiplex_res, allowed_mismatches)
  retained <- raw_filter_res$retained
  res <- raw_filter_res$demultiplex_res
  summary_res <- create_summary_res(retained,
                                    barcodes,
                                    res$assigned_barcodes,
                                    allowed_mismatches,
                                    mismatches)
  
  list(
    demultiplex_res = res,
    retained = retained,
    summary_res = summary_res
  )
}

validate_allowed_mismatches <- function(allowed_mismatches, barcodes) {
  n_barcode_sets <- length(barcodes)
  assert_that(
    length(allowed_mismatches) == 1L ||
      n_barcode_sets == length(allowed_mismatches),
    msg = "allowed_mismatches does not match the number of barcode sets"
  )
  if (length(allowed_mismatches) == 1L && n_barcode_sets > 1L) {
    allowed_mismatches <- rep(allowed_mismatches, n_barcode_sets)
  }
  names(allowed_mismatches) <- names(barcodes)
  allowed_mismatches
}

filter_sequences <- function(demultiplex_res, allowed_mismatches) {
  mismatches <- demultiplex_res$mismatches
  mismatches_above_threshold <- sweep(mismatches, 2L, allowed_mismatches, FUN = `>`)
  retained <- rowSums(mismatches_above_threshold) == 0L
  demultiplex_res$assigned_barcodes %<>% extract(retained, )
  demultiplex_res$mismatches %<>% extract(retained, )
  demultiplex_res$payload %<>% map(. %>% extract(retained))
  list(demultiplex_res = demultiplex_res, retained = retained)
}
