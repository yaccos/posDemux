#' Create a summary of match filtering
#' @description
#' \code{create_summary_res()} is a helper function in order to create a summary of the demultiplexing
#' and following match filtering. It is not designed to be invoked directly, but
#' its results will be returned automatically from
#' \code{\link{filter_demultiplex_res}}. This returned object has it own method
#' for printing the result in a user-friendly manner.
#' 
#' 
#' @inheritParams filter_demultiplex_res
#' @inheritParams combinatorial_demultiplex
#' @param barcodes A list of
#'  \code{\link[Biostrings:XStringSet-class]{XStringSet}} objects, the
#'  barcodes which were used for demultiplexing.
#' @param retained Logical vector with the same length as
#'  the number of reads in the input to the demultiplexer.
#'  \code{TRUE} if the corresponding read
#'  is retained. Corresponds to the field \code{retained} of the output of
#'  \code{\link{filter_demultiplex_res}}.
#' @param assigned_barcodes Character matrix of the assigned barcodes.
#' Corresponds to of the field \code{assigned_barcodes}
#' of \code{\link{combinatorial_demultiplex}}.
#' @param mismatches Integer matrix of the number of
#'  mismatches of each assigned barcode.
#' Corresponds to the field \code{mismatches} of
#'  \code{\link{combinatorial_demultiplex}}.
#'  
#'  
#' @return
#' \code{create_summary_res()} returns a list of S3 class \code{demultiplex_filter_summary}
#' providing diagnostics for the filtering process. It contains the
#' the following fields:
#' \itemize{
#' \item \code{n_reads}: The total number of reads in the dataset before filtering.
#' \item \code{n_removed}: The number of reads removed because demultiplexing failed.
#' \item \code{n_barcode_sets}: The number of barcode sets.
#' \item \code{n_barcode_combinations}: The possible number
#' of barcode combinations.
#' \item \code{n_unique_barcodes}: The number of observed unique barcode combinations
#' (i.e. features which may be cells) detected after filtering mismatches.
#' \item \code{n_estimated_features}: The estimated number of features having a 
#' detected combination of barcodes.
#' This number will always be greater or equal than \code{n_unique_barcodes} due
#' to barcode collisions.
#' \item \code{observed_collision_lambda}: The ratio of observed barcode
#'  combinations divided by the total number of possible barcode combinations.
#' \item \code{corrected_collision_lambda}: The ratio of estimated number of features
#' to the total number of possible barcode combinations.
#' \item \code{expected_collisions}: The statistically expected number
#' of barcode collisions or more precicely the expected number of
#' observed barcodes which correspond to two or more features.
#' \item \code{barcode_summary}: A list containing a summary for each barcode set.
#' Each element contains the following:
#' \itemize{
#' \item \code{width}: The width (number of nucleotides) of the barcode set.
#' \item \code{n_barcodes}: Number of query barcodes.
#' \item \code{n_allowed_mismatches}: Number of allowed mismatches for the barcode set.
#' \item \code{n_removed}: Number of reads having too many mismatches for this barcode set.
#' \item \code{mismatch_frame}: A \code{data.frame} with the two columns,
#' \code{n_mismatches} and \code{frequency} showing the number of reads for each
#' of the allowed number of mismatches for the given barcode set.
#' }
#' }
#' The \code{print()} method returns its output invisibly.
#' 
#' @details
#'  Following a uniform distribution of barcodes, the expected number
#'  of barcode collisions
#'  (observed barcodes combinations being composed of two or more features)
#'  is given by
#'  \deqn{N\left(1-e^{-\lambda}-\lambda e^{-\lambda}\right),}
#'  where \eqn{N} is the number of possible barcode combinations
#'  and \eqn{\lambda} is in this summary referred to as the collision lambda:
#'  \deqn{\lambda=\frac{n}{N},} where \eqn{n} is the number of features. 
#'  However, \eqn{n} is unknown as we cannot know how many features
#'  there were originally due to potential collisions.
#'  Utilizing the fact that the expected observed number of barcodes is given by
#'  \deqn{N\left(1-e^{-\lambda}\right),}
#'  we can correct the estimate for \eqn{\lambda} from the known value
#'  of the observed barcode combinations, and thus estimate the number
#'  of features and barcode collisions.
#'  
#'  While each unique feature can be conceptually thought of as single cell with its transcripts,
#'  realistic datasets have many features with relatively small numbers of reads which
#'  are artifacts and unlikely to correspond to true cells.
#'  
#' 
#' @example inst/examples/match_filter-examples.R
#' @export
create_summary_res <- function(retained,
                               barcodes,
                               assigned_barcodes,
                               allowed_mismatches,
                               mismatches) {
  allowed_mismatches <- validate_allowed_mismatches(allowed_mismatches, barcodes)
  n_removed <- sum(!retained)
  n_reads <- length(retained)
  n_barcodes_per_set <- map_int(barcodes, length)
  n_unique_barcodes <- assigned_barcodes %>%
    unique(MARGIN = 1L) %>%
    nrow()
  n_barcode_combinations <- prod(n_barcodes_per_set)
  n_estimated_features <- poisson_correct_n(n_barcode_combinations,
                                         n_unique_barcodes)
  n_barcodes_per_set <- map_int(barcodes, length)
  observed_collision_lambda <- n_unique_barcodes / n_barcode_combinations
  corrected_collision_lambda <- n_estimated_features / n_barcode_combinations
  expected_collisions <- poisson_estimate_collisions(n_barcode_combinations,
                                                     corrected_collision_lambda)
  barcode_summary <- imap(barcodes, function(barcode_set, barcode_name) {
    barcode_width <- width(barcode_set)[1L]
    n_allowed_mismatches <- allowed_mismatches[barcode_name]
    n_barcodes <- length(barcode_set)
    this_mismatch_vector <- mismatches[, barcode_name, drop =
                                         TRUE]
    
    mismatch_frame <- data.frame(n_mismatches =
                                   c(0L, seq_len(n_allowed_mismatches))) %>%
      mutate(frequency = outer(this_mismatch_vector, .data$n_mismatches, equals) %>%
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
    n_estimated_features = n_estimated_features,
    observed_collision_lambda = observed_collision_lambda,
    corrected_collision_lambda = corrected_collision_lambda,
    expected_collisions = expected_collisions,
    barcode_summary = barcode_summary
  )
  
  class(summary_res) <- "demultiplex_filter_summary"
  summary_res
}

poisson_estimate_collisions <- function(N, lambda) {
  N * (1 - exp(-lambda) - lambda * exp(-lambda))
}

# Estimates the true number of features taking barcode collisions into account
poisson_correct_n <- function(N, n_obs) {
  -N*log(1 - n_obs/ N)
}

#' @param x An object of class \code{demultiplex_filter_summary} from 
#' \code{create_summary_res()}.
#' @param ... Ignored
#' @importFrom purrr map_int walk2
#' @import glue
#' @rdname create_summary_res
#' @export
print.demultiplex_filter_summary <- function(x, ...) {
  glue("Total number of reads: {x$n_reads}") %>%
    cat("\n")
  removed_percentage <- x$n_removed / x$n_reads * 100
  glue(
    "Number of reads failing to demultiplex: \\
               {x$n_removed} ({removed_percentage %>% round(2L)}%)"
  ) %>%
    cat("\n")
  glue("Observed number of unique barcode combinations: {x$n_unique_barcodes}") %>%
    cat("\n")
  glue("Number of possible barcode combinations: {x$n_barcode_combinations}") %>%
    cat("\n")
  glue("Estimated number of features: {x$n_estimated_features %>% round(1L)}") %>%
    cat("\n")
  glue("Observed feature to barcode ratio: {x$observed_collision_lambda %>% signif(4L)}") %>% 
    cat("\n")
  glue("Corrected feature to barcode ratio: {x$corrected_collision_lambda %>% signif(4L)}") %>% 
    cat("\n")
  collision_percentage <- x$expected_collisions / x$n_unique_barcodes * 100
  glue(
    "Estimated number of observed barcode combinations
    corresponding to more than one feature: \\
    {x$expected_collisions %>% round(1L)} ({collision_percentage %>% round(2L)}%)"
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
