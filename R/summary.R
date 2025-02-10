#' Diagnostic demultiplexing results
#'
#' @param filter_results The output from \code{\link{filter_demultiplex_res}}
#' @param barcodes The list of \code{\link{XStringSet}} objects,
#' the barcodes which were used for demultiplexing.
#' @importFrom magrittr %>%
#' @import dplyr
#' @import glue
#' @returns NULL
#' @export
filter_summary <- function(filter_results, barcodes) {
  n_reads <- length(filter_results$retained)
  n_removed <- filter_results$n_removed
  n_barcodes_per_set <- map_int(barcodes, length)
  glue("Total number of reads: {filter_results$n_reads}") %>% 
    message()
  removed_percentage <- filter_results$n_removed/n_reads*100
  glue("Number of reads removed:\\
               {n_removed} ({removed_percentage %>% round(2L)}%)") %>% 
    message()
  n_barcode_sets <- length(filter_results$allowed_mismatches)
  glue("Number of barcode sets: {n_barcode_sets}") %>% message()
  barcode_sets <- colnames(filter_results$assigned_barcode)
  barcode_sets <- ifelse(is.null(barcode_sets), seq_len(n_barcode_sets), barcode_sets)
  for (barcode in barcode_sets) {
    cat(rep("-", 80L), "\n")
    glue("Barcode set: {barcode}") %>% message()
    this_removed <- filter_results$n_removed_per_barcode[i]
    this_percentage <- this_removed / n_reads * 100
    glue("Number of possible barcodes:\\
                 {n_barcodes_per_set[barcode]}") %>% 
      message()
    glue("Number of reads removed:\\
                 {this_removed} ({this_percentage %>% round(2L)}%)") %>% 
      message()
    n_allowed_mismatches <- filter_results$allowed_mismatches[barcode]
    for (mismatches in 0L:n_allowed_mismatches) {
      n_reads_with_mismatches <- filter_results$mismatches %>%
        extract(, barcode, drop=TRUE) %>%
        equals(mismatches) %>%
        sum()
      percentage <- n_reads_with_mismatches / n_reads * 100
      glue("Number of reads with {mismatches} mismatches:\\
                   {n_reads_with_mismatches} ({percentage %>% round(2L)}%)") %>% 
        message()
    }
  }
  cat(rep("-", 80L), "\n")
  n_unique_barcodes <- filter_results$assigned_barcodes %>%
    unique(MARGIN=1L) %>%
    nrow()
  glue("Number of unique barcodes detected: {n_unique_barcodes}") %>% 
    message()
  n_barcode_combinations <- prod(n_barcodes_per_set)
  glue("Number of possible barcode combinations: {n_barcode_combinations}") %>% 
    message()
  collision_lambda <- n_unique_barcodes / n_barcode_combinations
  expected_collisions <- n_unique_barcodes * collision_lambda
  collision_percentage <- collision_lambda * 100 
  glue("Expected number of barcode collisions:\\
       {expected_collisions}({collision_percentage %>% round(2L)}%)") %>% 
    message()
  return(NULL)
}

#'
#' Create a frequency table of the observed barcode combinations
#' 
#' @param filter_results The output from \code{\link{filter_demultiplex_res}} 
#' 
#'
#' @returns A data frame where each row correspond to a unique observed
#' barcode combination. The rows are sorted in descending order of frequency.
#' The first columns are specify the barcode assignment and the last columns were the
#' following:
#' \itemize{
#' \item \code{frequency}: The number of reads with the barcode combination
#' \item \code{cumulative_frequency}: The cumulative frequency of the barcode combination
#' \item \code{fractiion}: The fraction of reads with the barcode combination
#' \item \code{cumulative_fraction}: The cumulative fraction of the barcode combination
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#'
create_frequency_table <- function(filter_results) {
  filter_results$assigned_barcode %>%
    as.data.frame() %>%
    count(pick(everything()), name = "frequency", sort = TRUE) %>%
    mutate(cumulative_frequency = cumsum(.data$frequency),
           fraction = frequency / sum(.data$frequency) * 100,
           cumulative_fraction = cumulative_frequency / sum(.data$frequency))
}

frequency_plot <- function(frequency_table, sample_name ="",
                           binwidth = 1 / 20) {
  n_reads <- sum(frequency_table$frequency)
  ggplot(frequency_table, aes(x=.data$frequency)) +
    geom_histogram(binwidth = binwidth) +
    labs(title=glue("{sample_name}\\{ifelse(sample_name=='','',', ')}\\
                    Total reads = {n_reads}"),
         x="Frequency",
         y="log10(Number of reads)") +
    scale_y_log10()
}

knee_plot <- function(frequency_table, cutoff_line=NULL) {
    p <- ggplot(frequency_table, aes(x=n() %>% seq_len(),y=cumulative_fraction)) +
    labs(title=glue(""),
         x="Barcode (ordered largest to smallest)",
         y="Cumulative fraction of reads")
    if (!is.null(cutoff_line)) {
      p <- p + geom_vline(xintercept=cutoff_line, linetype="dashed")
    }
    p
}
