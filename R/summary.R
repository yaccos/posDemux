#' Diagnostic demultiplexing results
#'
#' @param filter_results The output from \code{\link{filter_demultiplex_res}}
#' @import glue
#' @returns NULL
#' @export
filter_summary <- function(filter_results) {
  n_reads <- length(filter_results$retained)
  message(glue("Total number of reads: {filter_results$n_reads}"))
  removed_percentage <- filter_results$n_removed/n_reads*100
  message(glue("Number of reads removed: {filter_results$n_removed} ({removed_percentage:.2f}%)"))
  n_barcode_sets <- length(filter_results$allowed_mismatches)
  message(glue("Number of barcode sets: {n_barcode_sets}"))
  barcode_sets <- colnames(filter_results$assigned_barcode)
  barcode_sets <- ifelse(is.null(barcode_sets), seq_len(n_barcode_sets), barcode_sets)
  for (i in seq_len(n_barcode_sets)) {
    this_removed <- filter_results$n_removed_per_barcode[i]
    this_percentage <- this_removed/n_reads*100
    message(glue("Number of reads removed for barcode set {i}: {this_removed} ({this_percentage:.2f}%)"))
  }
  
  
  
  
}