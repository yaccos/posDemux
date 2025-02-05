#' Filter demultiplexed reads based on the number of mismatches
#'
#' @param demultiplex_res Unprocessed output from
#'  \code{\link{combinatorial_demultiplex}}
#' @param allowed_mismatches Integer vector of length one or the same length
#' as the number of barcode segments. The maximum Hamming distance from including
#' a read in the output
#' @returns A list
#' @export
#'
filter_demultiplex_res <- function(demultiplex_res, allowed_mismatches){
  mismatches <- demultiplex_res$mismatches
  mismatches_above_threshold <- sweep(mismatches, 2L,
                                      allowed_mismatches, FUN = `>`)
  keep_sequence <- rowSums(mismatches_below_threshold) == 0L
  sequences_removed_per_barcode <- colSums(mismatches_above_threshold)
  
}
