#' Frequency table
#'
#' @description
#' Creates a sorted frequency table of each of the observed
#' barcode combinations. This function is indended to be used after running
#' [filter_demultiplex_res()] and before creating frequency plots,
#' knee plots, or selecting the number of barcodes to include.
#'
#' @param assigned_barcodes A character or integer matrix, corresponding to the field
#' \code{assigned_barcodes} from [combinatorial_demultiplex()]
#' or the field \code{demultiplex_res$assigned_barcodes}
#' from [filter_demultiplex_res()].
#'
#'
#' @returns A data frame where each row corresponds to a unique observed
#' barcode combination. The rows are sorted in descending order of frequency.
#' The first columns specify the
#' barcode assignment (e.g \code{bc3}, \code{bc2}, \code{bc1})
#' and the last columns were the following:
#' \itemize{
#' \item \code{frequency}: The number of reads with the barcode combination.
#' \item \code{cumulative_frequency}: The cumulative frequency of the barcode combination counted from the top.
#' \item \code{fraction}: The fraction of reads with the barcode combination.
#' \item \code{cumulative_fraction}: The cumulative fraction of the barcode combination counted from the top.
#' }
#' @import ggplot2
#' @importFrom rlang .data
#' @example inst/examples/freq_table-examples.R
#' @export
#'
create_freq_table <- function(assigned_barcodes) {
  assigned_barcodes %>%
    as.data.frame() %>%
    count(pick(everything()), name = "frequency", sort = TRUE) %>%
    mutate(
      cumulative_frequency = cumsum(.data$frequency),
      fraction = .data$frequency / sum(.data$frequency),
      cumulative_fraction = .data$cumulative_frequency / sum(.data$frequency)
    )
}

create_freq_table_from_encoding <- function(encoded_barcodes, mapping) {
  frequency <- table(encoded_barcodes) %>% sort(decreasing = TRUE)
  unique_encoded_barcodes <- names(frequency) %>% as.integer()
  unique_barcode_table <- decode(unique_encoded_barcodes, mapping) %>%
    mutate(
      frequency = frequency %>% as.vector() %>% unname(),
      cumulative_frequency = cumsum(.data$frequency),
      fraction = .data$frequency / sum(.data$frequency),
      cumulative_fraction = .data$cumulative_frequency / sum(.data$frequency)
    )
}

create_freq_table_from_count_table <- function(count_table, mapping) {
  # Materializes the Rcpp count table
  count_res <- get_count_table(count_table)
  decode(count_res$encoding, mapping) %>%
    mutate(frequency = count_res$frequency) %>%
    arrange(desc(.data$frequency)) %>%
    mutate(
      cumulative_frequency = cumsum(.data$frequency),
      fraction = .data$frequency / sum(.data$frequency),
      cumulative_fraction = .data$cumulative_frequency / sum(.data$frequency)
    )
}

#' Diagnostic plots from demultiplexing
#'
#' @description
#' Diagnostic plots for determining the effect of the barcode cutoff.
#' \code{freq_plot()} shows a histogram or distribution plot of the
#' number of reads for each barcode
#' combination, whereas \code{knee_plot()} shows the cumulative fraction of reads
#' ranked by the frequency of the barcode combinations in descending order.
#'
#' @param freq_table The frequency table
#'  from [create_freq_table()].
#'
#' @param cutoff Optional scalar numeric, the
#' x-coordinate for drawing a vertical dashed line
#' in the plots in order to indicate the cutoff. Please note that this argument
#' is interpreted literally, meaning that in order to correctly display
#' the same cutoff on both type of plots,
#' the cutoff value has to be transformed. In order to safely convert between
#' the two types of cutoffs, use the functions [bc_to_freq_cutoff()]
#' and [freq_to_bc_cutoff()].
#' @param type The type of frequency plot to make, either \code{"histogram"}
#' or \code{"density"}.
#' @param log_scale_x Logical: Should a log scale be applied to the x-axis of the
#' frequency plot?
#' @param log_scale_y Logical: Should a log scale be applied to the y-axis of the
#' frequency plot?
#' @param scale_by_reads Logical: Should the y-axis of the plot be scaled by the
#' number of reads on the x-axis?
#'
#' @seealso [bc_to_freq_cutoff()] [freq_to_bc_cutoff()]
#'
#' @importFrom rlang .data
#'
#' @returns A \code{\link[ggplot2]{ggplot}} object which can be displayed immediately or
#' further modified.
#'
#' @example inst/examples/freq_plot-examples.R
#' @export
freq_plot <- function(freq_table,
                      cutoff = NULL,
                      type = "histogram",
                      log_scale_x = TRUE,
                      log_scale_y = FALSE,
                      scale_by_reads = FALSE) {
  n_reads <- sum(freq_table$frequency)
  if (!scale_by_reads) {
    plot_type <- switch(type,
      histogram = \() geom_histogram(),
      density = \() stat_density()
    )
  } else {
    plot_type <- switch(type,
      histogram = \() geom_histogram(aes(y = after_stat(.data$count * .data$x))),
      density = \() stat_density(aes(y = after_stat(.data$count * .data$x / sum(.data$count))))
    )
  }

  if (is.null(plot_type)) {
    stop("The type argument must either be 'histogram' or 'density'")
  }
  if (type == "histogram") {
    ylab <- "Frequency"
  } else if (log_scale_y) {
    # This plot is special since the values closest to zero
    # represent the most reads
    ylab <- "log10(Relative Frequency)"
  } else {
    ylab <- "Relative Frequency"
  }
  p <- ggplot(freq_table, aes(x = .data$frequency)) +
    plot_type() +
    labs(x = "Number of reads", y = ylab)
  if (log_scale_x) {
    p <- p + scale_x_log10()
  }
  if (log_scale_y) {
    p <- p + scale_y_log10()
  }
  if (!is.null(cutoff)) {
    p <- p + geom_vline(xintercept = cutoff, linetype = "dashed")
  }
  p
}

#' @rdname freq_plot
#' @import dplyr
#' @importFrom rlang .data
#' @export
knee_plot <- function(freq_table, cutoff = NULL) {
  augmented_frequency_table <- freq_table %>%
    select("cumulative_fraction") %>%
    mutate(index = seq_len(n())) %>%
    add_row(index = 0L, cumulative_fraction = 0)
  p <- ggplot(
    augmented_frequency_table,
    aes(x = .data$index, y = .data$cumulative_fraction)
  ) +
    geom_line() +
    labs(x = "Barcode (ordered largest to smallest)", y = "Cumulative fraction of reads") +
    ylim(0, 1) +
    xlim(0L, nrow(augmented_frequency_table))
  if (!is.null(cutoff)) {
    p <- p + geom_vline(xintercept = cutoff, linetype = "dashed")
  }
  p
}
