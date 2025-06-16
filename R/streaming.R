#' Demultiplexing with streaming
#'
#' @description
#' This function provides an interface to
#' \code{\link{combinatorial_demultiplex}} and
#' \code{\link{filter_demultiplex_res}} such that reads are streamed in chunks
#' instead having to load everything at once, hence reducing memory consumption.
#' It accepts two functions which are called once per chunk.
#' A data loader function for producing the sequences of the chunk and an
#' archiver writing the results to file. This function also decides the size of the chunk.
#' for sequences and an online archiver function
#'
#' @param state_init The initial state to pass into \code{loader}
#' @param loader Function loading the reads. It has the signature
#' \code{f(state)}, where \code{state} is a user-defined object
#' which is initialized to be \code{state_init} and for the subsequent
#' iterations taken as the \code{state} field of the output of \code{archiver}.
#' Its return value is a list with the following fields:
#' \itemize{
#' \item \code{state}: The state to be passed into \code{archiver}
#' \item \code{sequences}: A \code{\link{XStringSet}} object, the sequences
#' to be demultiplexed in the current chunk.
#' \item \code{should_terminate}: A scalar logical. If \code{TRUE}, the demultiplexing
#' process terminates and the final results are returned.
#' Notice that this termination happens before the sequences of the final
#' call to \code{loader} are demultiplexed
#' }
#'
#' @param archiver Function taking care of archiving the demultiplexing results.
#' Its arguments are:
#' \itemize{
#' \item \code{state}: The state of the process returned by \code{loader}
#' \item \code{filtered_res}: The output from running
#' \code{\link{combinatorial_demultiplex}} and
#' \code{\link{filter_demultiplex_res}} on the data expect that the
#' field \code{summary_res} is missing.
#' }
#' Its output is a state object fed into the next call to \code{loader}.
#' @details
#' The data loader decides the size of each chunk.
#' While this framework does not provide any restriction on the \code{state}
#' object, the loader and archiver must be written such that the state objects
#' they return are compatible.
#' Since the data loader alone decides when to terminate,
#' bad terminations crieria can cause a runaway loop.
#'
#' @seealso [filter_demultiplex_res()] and [combinatorial_demultiplex()]
#' @inheritParams combinatorial_demultiplex
streaming_demultiplex <- function(state_init,
                                  loader,
                                  archiver,
                                  barcodes,
                                  allowed_mismatches,
                                  segments,
                                  segment_lengths) {
  summary <- summary_init(barcodes, allowed_mismatches)
  loader_res <- loader(state_init)
  encoded_barcodes <- create_encoding_vector()
  barcode_table <- map(barcodes, names)
  barcode_mapping <- get_mapping(barcode_table)
  while (!loader_res$should_terminate) {
    state <- loader_res$state
    sequences <- loader_res$sequences
    demultiplex_res <- combinatorial_demultiplex(sequences, barcodes, segments, segment_lengths)
    filtered_res <- filter_sequences(demultiplex_res, allowed_mismatches)
    summary <- summary_update(summary, filtered_res$retained,
                              filtered_res$demultiplex_res$mismatches)
    this_barcode_encoding <- as.data.frame(filtered_res$demultiplex_res$assigned_barcodes,
                                           row.names = FALSE) %>% encode(barcode_mapping)
    grow_encoding_vector(encoded_barcodes, this_barcode_encoding)
    state <- archiver(state, filtered_res)
    loader_res <- loader(state)
  }
  final_assigned_barcodes <- get_encoding_vector(encoded_barcodes,
                                                 barcode_mapping) %>% 
    decode(encoded_barcodes, barcode_mapping)
  freq_table <- final_assigned_barcodes %>%
    as.matrix() %>% 
    create_frequency_table()
  final_summary <- summary_finalize(summary)
  list(filter_summary=final_summary, freq_table=freq_table)
}

summary_init <- function(barcodes, allowed_mismatches) {
  n_removed <- 0L
  n_reads <- 0L
  n_barcodes_per_set <- map_int(barcodes, length)
  n_barcode_combinations <- prod(n_barcodes_per_set)
  n_barcodes_per_set <- map_int(barcodes, length)
  barcode_summary <- imap(barcodes, function(barcode_set, barcode_name) {
    barcode_width <- width(barcode_set)[1L]
    n_allowed_mismatches <- allowed_mismatches[barcode_name]
    n_barcodes <- length(barcode_set)
    this_mismatch_vector <- mismatches[, barcode_name, drop =
                                         TRUE]
    mismatch_frame <- data.frame(n_mismatches =
                                   c(0L, seq_len(n_allowed_mismatches)),
                                 frequency = 0L)
    list(
      width = barcode_width,
      n_barcodes = n_barcodes,
      n_allowed_mismatches = n_allowed_mismatches,
      n_removed = 0L,
      mismatch_frame = mismatch_frame
    )
  })
  
  list(
    n_reads = n_reads,
    n_removed = n_removed,
    n_barcode_sets = barcodes %>% length(),
    n_barcode_combinations = n_barcode_combinations,
    barcode_summary = barcode_summary
  )
}


summary_update <- function(filter_summary, retained_sequences, mismatches) {
  within(filter_summary, {
    n_reads <- n_reads + length(retained_sequences)
    n_removed <- n_removed + sum(!retained_sequences)
    barcode_summary <- imap(barcode_summary, function(barcode_summary, barcode_name) {
      this_mismatch_vector <- mismatches[, barcode_name, drop = TRUE]
      within(barcode_summary, {
        n_removed <- n_removed + sum(this_mismatch_vector > n_allowed_mismatches)
        mismatch_frame <- within(mismatch_frame, {
          # Yes, this is nested 5 levels deep, but currently I not know of any
          # solution more aestetic
          frequency <- frequency +
            outer(this_mismatch_vector, n_mismatches, equals) %>%
            colSums()
        })
      })
    })
  })
}

summary_finalize <- function(filter_summary, n_unique_barcodes) {
  filter_summary <- within(filter_summary, {
    n_unique_barcodes <- n_unique_barcodes
    n_estimated_cells <- poisson_correct_n(n_barcode_combinations, n_unique_barcodes)
    observed_collision_lambda <- n_unique_barcodes / n_barcode_combinations
    corrected_collision_lambda <- n_estimated_cells / n_barcode_combinations
    expected_collisions <- poisson_estimate_collisions(n_barcode_combinations, corrected_collision_lambda)
  })
  class(filter_summary) <- "demultiplex_filter_summary"
  filter_summary
}
