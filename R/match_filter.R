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
#' \item \code{summary_res}: A list of S3 class \code{demultiplex_filter_summary} 
#' providing diagnostics for the filtering process. It contains the
#' the following fields:
#' \itemize{
#' \item \code{n_reads}: The total number of reads in the dataset before filtering.
#' \item \code{n_removed}: The number of reads removed by filtering.
#' \item \code{n_barcode_sets}: The number of barcode sets.
#' \item \code{n_barcode_combinations}: The possible number
#' of barcode combinations.
#' 
#' \item \code{n_unique_barcodes}: The number of unique barcode combinations
#' (i.e. cells) detected after filtering.
#' \item \code{collision_lambda}: The expected relative frequency of barcode collisions.
#' \item \code{expected_collisions}: The statistically expected number of barcode collisions.
#' \item \code{barcode_summary}: A list containing a summary for each barcode set.
#' Each element contains the following:
#' \itemize{
#' \item \code{width}: The width (number of nucleotide) of the barcode set.
#' \item \code{n_barcodes}: Number of query barcodes.
#' \item \code{n_allowed_mismatches}: Number of allowed mismatches for the barcode set.
#' \item \code{n_removed}: Number of reads having too many mismatches for this barcode set.
#' \item \code{mismatch_frame}: A \code{data.frame} with the two columns, 
#' \item \code{n_mismatches} and \code{frequency} showing the number of reads for each
#' of the allowed number of mismatches for the given barcode set.
#' }
#' }
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
  assert_that(length(allowed_mismatches) == 1L || 
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
  # n_sequences_removed_per_barcode <- colSums(mismatches_above_threshold)
  n_removed <- sum(!retained_sequences)
  n_reads <- length(retained_sequences)
  n_barcodes_per_set <- map_int(barcodes, length)
  n_unique_barcodes <- res$assigned_barcodes %>%
    unique(MARGIN=1L) %>%
    nrow()
  n_barcode_combinations <- prod(n_barcodes_per_set)
  n_barcodes_per_set <- map_int(barcodes, length)
  collision_lambda <- n_unique_barcodes / n_barcode_combinations
  expected_collisions <- n_unique_barcodes * collision_lambda
  barcode_summary <- imap(demultiplex_res$barcodes,
                          function(barcode_set, barcode_name) {
                            barcode_width <- width(barcode_set)[1L]
                            n_allowed_mismatches <- allowed_mismatches[barcode_name]
                            n_barcodes <- length(barcode_set)
                            this_mismatch_vector <- mismatches[, barcode_name, drop=TRUE]
                            
                            mismatch_frame <- data.frame(n_mismatches = 
                                                           c(0L, seq_len(n_allowed_mismatches))
                            ) %>%
                              mutate(frequency = outer(this_mismatch_vector, n_mismatches, equals) %>% 
                                       colSums()
                              )
                            this_removed <-  sum(this_mismatch_vector > n_allowed_mismatches)
                            list(
                                 width = barcode_width,
                                 n_barcodes = n_barcodes,
                                 n_allowed_mismatches = n_allowed_mismatches,
                                 n_removed = this_removed,
                                 mismatch_frame = mismatch_frame
                                 )
                          }
  )
  
  summary_res <- list(
    n_reads=n_reads,n_removed = n_removed,
    n_barcode_sets=n_barcode_sets,
    n_barcode_combinations = n_barcode_combinations,
    n_unique_barcodes = n_unique_barcodes,
    collision_lambda = collision_lambda,
    expected_collisions = expected_collisions,
    barcode_summary = barcode_summary
  )
  
  
  class(summary_res) <- "demultiplex_filter_summary"
  
  list(demultiplex_res = res, retained = retained_sequences,
       summary_res = summary_res)
}

filter_sequences <- function(demultiplex_res, allowed_mismatches) {
  mismatches <- demultiplex_res$mismatches
  mismatches_above_threshold <- sweep(mismatches, 2L,
                                      allowed_mismatches, FUN = `>`)
  retained_sequences <- rowSums(mismatches_above_threshold) == 0L
  demultiplex_res$assigned_barcodes %<>% extract(retained_sequences,)
  demultiplex_res$mismatches %<>% extract(retained_sequences,)
  demultiplex_res$payload %<>% extract(retained_sequences)
  list(res=demultiplex_res, retained_sequences=retained_sequences)
}

#' Diagnostic demultiplexing results
#' 
#' @description
#' Prints diagnostic information about the results of demultiplexing
#'  and subsequent filtering, including the results per barcode set. 
#' 
#'
#' @param x The field \code{filter_summary}  from \code{\link{filter_demultiplex_res}}
#' @param ... Ignored
#' @importFrom magrittr %>% extract equals
#' @importFrom purrr map_int
#' @import dplyr
#' @import glue
#' @returns Its input, invisibly.
#' @export
print.demultiplex_filter_summary <- function(x, ...) {
  glue("Total number of reads: {x$n_reads}") %>% 
    cat("\n")
  removed_percentage <- x$n_removed/x$n_reads * 100
  glue("Number of reads removed by filtering: \\
               {x$n_removed} ({removed_percentage %>% round(2L)}%)") %>% 
    cat("\n")
  glue("Number of unique barcodes detected: {x$n_unique_barcodes}") %>% 
    cat("\n")
  glue("Number of possible barcode combinations: {x$n_barcode_combinations}") %>% 
    cat("\n")
  collision_percentage <- x$collision_lambda * 100 
  glue("Expected number of barcode collisions: \\
       {x$expected_collisions} ({collision_percentage %>% round(2L)}%)") %>% 
    cat("\n")
  glue("Number of barcode sets: {x$n_barcode_sets}") %>% cat("\n")
  iwalk(x$barcode_summary, function(res,barcode_name) {
    cat(rep("-", 80L), "\n")
    glue("Barcode set: {barcode_name}") %>% cat("\n")
    glue("Barcode width: {res$width}") %>% cat("\n")
    glue("Number of possible barcodes: \\
                 {res$n_barcodes}") %>% 
      cat("\n")
    glue("Number of allowed mismatches: {res$n_allowed_mismatches}") %>% cat("\n")
    walk2(res$mismatch_frame$n_mismatches, res$mismatch_frame$frequency, 
         function(mismatches, frequency) {
           percentage <- frequency / x$n_reads * 100
           glue("Number of reads with {mismatches} mismatches: \\
                   {frequency} ({percentage %>% round(2L)}%)") %>% 
             cat("\n")
         }
           )
    percentage <- res$n_removed / x$n_reads * 100
    glue("Number of reads above mismatch threshold: \\
                 {res$n_removed} ({percentage %>% round(2L)}%)") %>% 
      cat("\n")
  }
  )
  cat(rep("-", 80L), "\n")
  invisible(x)
}
