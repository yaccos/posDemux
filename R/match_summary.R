create_summary_res <- function(retained_sequences,
                               barcodes,
                               assigned_barcodes,
                               allowed_mismatches,
                               mismatches) {
  n_removed <- sum(!retained_sequences)
  n_reads <- length(retained_sequences)
  n_barcodes_per_set <- map_int(barcodes, length)
  n_unique_barcodes <- assigned_barcodes %>%
    unique(MARGIN = 1L) %>%
    nrow()
  n_barcode_combinations <- prod(n_barcodes_per_set)
  n_barcodes_per_set <- map_int(barcodes, length)
  collision_lambda <- n_unique_barcodes / n_barcode_combinations
  expected_collisions <- n_unique_barcodes * collision_lambda / 2
  barcode_summary <- imap(barcodes, function(barcode_set, barcode_name) {
    barcode_width <- width(barcode_set)[1L]
    n_allowed_mismatches <- allowed_mismatches[barcode_name]
    n_barcodes <- length(barcode_set)
    this_mismatch_vector <- mismatches[, barcode_name, drop =
                                         TRUE]
    
    mismatch_frame <- data.frame(n_mismatches =
                                   c(0L, seq_len(n_allowed_mismatches))) %>%
      mutate(frequency = outer(this_mismatch_vector, n_mismatches, equals) %>%
               colSums())
    this_removed <-  sum(this_mismatch_vector > n_allowed_mismatches)
    list(
      width = barcode_width,
      n_barcodes = n_barcodes,
      n_allowed_mismatches = n_allowed_mismatches,
      n_removed = this_removed,
      mismatch_frame = mismatch_frame
    )
  })
  
  summary_res <- list(
    n_reads = n_reads,
    n_removed = n_removed,
    n_barcode_sets = barcodes %>% length(),
    n_barcode_combinations = n_barcode_combinations,
    n_unique_barcodes = n_unique_barcodes,
    collision_lambda = collision_lambda,
    expected_collisions = expected_collisions,
    barcode_summary = barcode_summary
  )
  
  class(summary_res) <- "demultiplex_filter_summary"
  summary_res
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
  removed_percentage <- x$n_removed / x$n_reads * 100
  glue(
    "Number of reads removed by filtering: \\
               {x$n_removed} ({removed_percentage %>% round(2L)}%)"
  ) %>%
    cat("\n")
  glue("Number of unique barcodes detected: {x$n_unique_barcodes}") %>%
    cat("\n")
  glue("Number of possible barcode combinations: {x$n_barcode_combinations}") %>%
    cat("\n")
  collision_percentage <- x$collision_lambda * 100
  glue(
    "Expected number of barcode collisions: \\
       {x$expected_collisions} ({collision_percentage %>% round(2L)}%)"
  ) %>%
    cat("\n")
  glue("Number of barcode sets: {x$n_barcode_sets}") %>% cat("\n")
  iwalk(x$barcode_summary, function(res, barcode_name) {
    cat(rep("-", 80L), "\n")
    glue("Barcode set: {barcode_name}") %>% cat("\n")
    glue("Barcode width: {res$width}") %>% cat("\n")
    glue("Number of possible barcodes: \\
                 {res$n_barcodes}") %>%
      cat("\n")
    glue("Number of allowed mismatches: {res$n_allowed_mismatches}") %>% cat("\n")
    walk2(res$mismatch_frame$n_mismatches, res$mismatch_frame$frequency, function(mismatches, frequency) {
      percentage <- frequency / x$n_reads * 100
      glue(
        "Number of reads with {mismatches} mismatches: \\
                   {frequency} ({percentage %>% round(2L)}%)"
      ) %>%
        cat("\n")
    })
    percentage <- res$n_removed / x$n_reads * 100
    glue(
      "Number of reads above mismatch threshold: \\
                 {res$n_removed} ({percentage %>% round(2L)}%)"
    ) %>%
      cat("\n")
  })
  cat(rep("-", 80L), "\n")
  invisible(x)
}
