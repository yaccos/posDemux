#' Filter demultiplexed reads based on the number of mismatches
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
#' as the number of barcode segments. The maximum Hamming distance from including
#' a read in the output
#' @returns A list with the following elements:
#' \itemize{
#' \item \code{demultiplex_res}: The contents of \code{demultiplex_res} with
#' the sequences filtered
#' \item \code{n_removed}: Integer, the total number of sequences removed
#' \item \code{n_removed_per_barcode}: Integer vector, 
#' the number of sequences removed per barcode
#' 
#' \item \code{retained}: Logical vector with the same length as
#'  the number of reads in the input. \code{TRUE} if the corresponding read
#'  is retained.
#'  Useful for future filtering of paired-end reads.
#'  \item \code{allowed_mismatches}: Integer vector, the number
#'  of allowed mismatches for each barcode set
#'  }
#' @details
#' The value of \code{n_removed} does not in general equal the sum of of
#' \code{n_removed_per_barcode} since a read can have too many mismatches
#' with multiple barcodes. 
#' 
#' @export
#'
filter_demultiplex_res <- function(demultiplex_res, allowed_mismatches) {
  mismatches <- demultiplex_res$mismatches
  mismatches_above_threshold <- sweep(mismatches, 2L,
                                      allowed_mismatches, FUN = `>`)
  retained_sequences <- rowSums(mismatches_below_threshold) == 0L
  n_sequences_removed_per_barcode <- colSums(mismatches_above_threshold)
  n_removed_sequences <- sum(!retained_sequences)
  res <- demultiplex_res
  res$assigned_barcodes <- res$assigned_barcodes[retained_sequences,]
  res$mismatches <- res$mismatches[retained_sequences,]
  res$payload <- res$payload[retained_sequences]
  list(demultiplex_res=res, n_removed=n_removed_sequences,
       n_removed_per_barcode=n_sequences_removed_per_barcode,
       retained=retained_sequences,allowed_mismatches)
}
