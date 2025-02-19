#' Diagnostic demultiplexing results
#'
#' @param filter_results The output from \code{\link{filter_demultiplex_res}}
#' @param barcodes The list of \code{\link{XStringSet}} objects,
#' the barcodes which were used for demultiplexing.
#' @importFrom magrittr %>% extract equals
#' @importFrom purrr map_int
#' @import dplyr
#' @import glue
#' @returns NULL
#' @export
filter_summary <- function(filter_results, barcodes) {
  n_reads <- length(filter_results$retained) + filter_results$n_removed
  n_removed <- filter_results$n_removed
  n_barcodes_per_set <- map_int(barcodes, length)
  glue("Total number of reads: {n_reads}") %>% 
    message()
  removed_percentage <- filter_results$n_removed/n_reads*100
  glue("Number of reads removed: \\
               {n_removed} ({removed_percentage %>% round(2L)}%)") %>% 
    message()
  n_barcode_sets <- length(filter_results$allowed_mismatches)
  glue("Number of barcode sets: {n_barcode_sets}") %>% message()
  barcode_sets <- colnames(filter_results$demultiplex_res$assigned_barcode)
  if (is.null(barcode_sets)) {
    barcode_sets <- seq_len(n_barcode_sets)
  }
  for (barcode in barcode_sets) {
    cat(rep("-", 80L), "\n")
    glue("Barcode set: {barcode}") %>% message()
    this_removed <- filter_results$n_removed_per_barcode[barcode]
    this_percentage <- this_removed / n_reads * 100
    glue("Number of possible barcodes: \\
                 {n_barcodes_per_set[barcode]}") %>% 
      message()
    glue("Number of reads removed: \\
                 {this_removed} ({this_percentage %>% round(2L)}%)") %>% 
      message()
    n_allowed_mismatches <- filter_results$allowed_mismatches[barcode]
    for (mismatches in 0L:n_allowed_mismatches) {
      n_reads_with_mismatches <- filter_results$demultiplex_res$mismatches %>%
        extract(, barcode, drop=TRUE) %>%
        equals(mismatches) %>%
        sum()
      percentage <- n_reads_with_mismatches / n_reads * 100
      glue("Number of reads with {mismatches} mismatches: \\
                   {n_reads_with_mismatches} ({percentage %>% round(2L)}%)") %>% 
        message()
    }
  }
  cat(rep("-", 80L), "\n")
  n_unique_barcodes <- filter_results$demultiplex_res$assigned_barcodes %>%
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
  glue("Expected number of barcode collisions: \\
       {expected_collisions} ({collision_percentage %>% round(2L)}%)") %>% 
    message()
  invisible(NULL)
}

#'
#' Frequency table 
#' 
#' @description
#' Creates a sorted frequency table of each of the observed
#' barcode combinations. This function is indended to be used after running
#' \code{\link{filter_demultiplex_res}} and before creating frequency plots,
#' knee plots or selecting the number of barcodes to include.
#' 
#' @param filter_results A character or integer matrix, corresponding to the field
#' \code{assigned_barcode} from \code{\link{filter_demultiplex_res}} 
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
#' }
#' @import ggplot2
#' @importFrom rlang .data
#' @export
#'
create_frequency_table <- function(filter_results) {
  filter_results %>%
    as.data.frame() %>%
    count(pick(everything()), name = "frequency", sort = TRUE) %>%
    mutate(cumulative_frequency = cumsum(.data$frequency),
           fraction = frequency / sum(.data$frequency),
           cumulative_fraction = cumulative_frequency / sum(.data$frequency))
}

#' Diagnostic plots from demultiplexing
#' 
#' @description
#' Diagnostic plots for determining the effect of the barcode cutoff.
#' \code{frequency_plot()} shows a histogram of the number of reads for each barcode
#' combination, whereas \code{knee_plot()} shows the cumulative fraction of reads
#' ranked by the frequency of the barcode combinations in descending order.
#' 
#' @param frequency_table The frequency table
#'  from \code{\link{create_frequency_table}}
#'  
#' @param width Optional positive scalar numeric,
#' the width of of the bars in the frequency plot
#' @param cutoff Optional scalar numeric, the
#' x-coordinate for drawing a vertical dashed line
#' in the plots in order to indicate the cutoff. Please note that this argument
#' is interpreted literally, meaning the in order to correctly display
#' the same cutoff on both type of plots,
#' the cutoff values has to be transformed. In order to safely convert between
#' the two types of cutoffs, use the functions
#' \code{\link{bc_to_frequency_cutoff}}
#' and \code{\link{frequency_to_bc_cutoff}}.
#' 
#' @seealso [bc_to_frequency_cutoff()] [frequency_to_bc_cutoff()]
#'  
#' @returns A \code{\link{ggplot}} which can be displayed immediately or
#' further modified
#'
#' @export
frequency_plot <- function(frequency_table,
                           width = NULL, cutoff = NULL) {
  n_reads <- sum(frequency_table$frequency)
  p <- ggplot(frequency_table, aes(x=.data$frequency)) +
    stat_count() +
    labs(x="Number of reads", y="Frequency") +
    xlim(0, max(frequency_table$frequency) + 1)
    # scale_y_log10()
  if(!is.null(cutoff)) {
    p <- p + geom_vline(xintercept=cutoff, linetype="dashed")
  }
  p
}

#' @rdname frequency_plot
#' @export
knee_plot <- function(frequency_table, cutoff=NULL) {
    augmented_frequency_table <- frequency_table %>%
      select(cumulative_fraction) %>%
      mutate(index = seq_len(n())) %>% 
      add_row(index=0L, cumulative_fraction=0)
    p <- ggplot(augmented_frequency_table,
                aes(x=index, y=cumulative_fraction)) +
      geom_line() +
    labs(
      x="Barcode (ordered largest to smallest)",
      y="Cumulative fraction of reads"
      ) +
      ylim(0,1) +
      xlim(0L, nrow(augmented_frequency_table))
    if (!is.null(cutoff)) {
      p <- p + geom_vline(xintercept=cutoff, linetype="dashed")
    }
    p
}
