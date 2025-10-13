log_progress <- function(msg) {
  message(glue("{date()} => {msg}"))
}

#' Suggested setup for FASTQ streaming
#' @description
#' Even though the user can define the arguments \code{state_init},
#' \code{loader}, and \code{archiver}
#' for [streaming_demultiplex()], this approach is only recommended
#' for advanced users. This functions defines a premade combinations of
#' these three arguments which should be suitable in most cases.
#' The loader streams a FASTQ file in chunks using [ShortRead::FastqStreamer()]
#' and the archiver outputs
#' a data frame to file consisting of the following of the read name (\code{read}),
#' the sequences of all payloads (e.g. \code{UMI}), and barcode assignments
#' (\code{c("bc3","bc2","bc1")}).
#'
#' @param input_file The path to the FASTQ file to be used for demultiplexing
#' @param output_table_file The path to which the output
#' barcode table will be written
#' @param chunk_size Integer, the number of reads to process in each chunk
#' @param verbose Logical scalar: Should the progress be displayed?
#' @param min_width Optional integer scalar: Minimum width of the sequences to keep.
#' For reads which are shorter than this, a warning it emitted and the
#' reads are removed and ignored and thus not appear in any statistics.
#' The data loader is **not** supposed to be used as a length filter, so this option
#' is more like an escape hatch for being able to deal with sequences which have not been
#' properly filtered beforehand.
#'  
#' 
#' @details
#' If the read names have any spaces in them, 
#' the loader will only keep the portion of
#' the read name preceding the first space. This is due to the Illumina
#' platform's behavior of encoding the sequencing direction (forward or reverse)
#' past the space. 
#' Keeping the read names with the space is usually not desirable as it makes
#' the resulting barcode table more confusing and makes it more difficult to
#' group the forward and reverse reads together afterwards.
#' 
#'
#' @returns A list with the following elements, all of which are intended to be
#' used as the corresponding arguments to [streaming_demultiplex()]:
#' \itemize{
#' \item \code{state_init}
#' \item \code{loader}
#' \item \code{archiver}
#' }
#' @importFrom purrr list_cbind
#' @export
streaming_callbacks <- function(input_file,
                                output_table_file,
                                chunk_size = 1e6,
                                verbose = TRUE,
                                min_width = NULL) {
  res <- list()
  res$state_init <- list(
    total_reads = 0L,
    demultiplexed_reads = 0L,
    output_table_initialized = FALSE
  )
  
  res$loader <- function(state) {
    if (!state$output_table_initialized) {
      if(verbose) log_progress("Initializing FASTQ stream and output table")
      state$istream <- ShortRead::FastqStreamer(input_file, n = chunk_size)
    }
    raw_chunk <- ShortRead::yield(state$istream)
    chunk  <- ShortRead::sread(raw_chunk)
    # For pair-end reads, we usually don't want what is trailing after the space
    # in order to have the same identifiers to both forward and reverse reads 
    names(chunk) <- ShortRead::id(raw_chunk) %>% {sub(" .*$", "", .)}
    n_reads_in_chunk <- length(chunk)
    if (n_reads_in_chunk == 0L && state$output_table_initialized) {
      # The case when the initial chunk is empty is given special treatment since
      # we want to create the table regardless
      # No more reads to process, make the outer framework return
      #Clean up resources
      close(state$istream)
      state$istream <- NULL
      final_res <- list(
        state = state,
        sequences = NULL,
        should_terminate = TRUE
      )
      if(verbose) log_progress("Done demultiplexing")
      return(final_res)
    }
    if (min_width %>% is.null() %>% magrittr::not()){
      chunk <- warn_sufficient_length(chunk, min_width)
      n_reads_in_chunk <- length(chunk)
    }
    state$total_reads <- state$total_reads + n_reads_in_chunk
    list(
      state = state,
      sequences = chunk,
      should_terminate = FALSE
    )
  }
  
  res$archiver <- function(state, filtered_res) {
    barcode_matrix  <- filtered_res$demultiplex_res$assigned_barcodes
    barcode_names <- colnames(barcode_matrix)
    read_names <- rownames(barcode_matrix)
    # If the table has no rows, we may risk getting a NULL value
    if (is.null(read_names)) {
      read_names <- character()
    }
    barcode_table <- as.data.frame(barcode_matrix)
    read_name_table <- data.frame(read = read_names)
    payload_table <- purrr::map(filtered_res$demultiplex_res$payload, as.character) %>%
      {
        rlang::exec(data.frame, !!!.)
      }
    
    chunk_table <- cbind(read_name_table, payload_table, barcode_table)
    if (!state$output_table_initialized) {
      append <- FALSE
      state$output_table_initialized <- TRUE
    } else {
      append <- TRUE
    }
    data.table::fwrite(
      x = chunk_table,
      file = output_table_file,
      append = append,
      row.names = FALSE,
      col.names = !append,
      sep = "\t",
      eol = "\n"
    )
    state <- within(state, {
      demultiplexed_reads <- demultiplexed_reads + nrow(barcode_matrix)
      if(verbose)
        log_progress(
          "Processed {total_reads} reads, successfully demultiplexed {demultiplexed_reads} reads so far..."
        )
    })
    state
  }
  res
}
