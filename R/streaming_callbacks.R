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
#' @param chunk_size The number of reads to process in each chunk
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
streaming_callbacks <- function(input_file, output_table_file,
                          chunk_size=1e6) {
  res <- list()
  res$state_init <- list(total_reads=0L, demultiplexed_reads=0L,
                     output_table_initialized=FALSE)
  
  res$loader <- function(state) {
    if (!state$output_table_initialized) {
      message("Initializing FASTQ stream and output table")
      state$istream <- ShortRead::FastqStreamer(input_file,
                                                        n = chunk_size)
    }
    chunk  <- ShortRead::yield(state$istream) %>% ShortRead::sread()
    n_reads_in_chunk <- length(chunk)
    if (n_reads_in_chunk == 0L && state$output_table_initialized) {
      # The case when the initial chunk is empty is given special treatment since
      # we want to create the table regardless
      # No more reads to process, make the outer framework return
      #Clean up resources
      close(state$istream)
      state$istream <- NULL
      final_res <- list(state = state, sequences=NULL, should_terminate=TRUE)
      message(glue("Done demultiplexing {state$total_reads} reads"))
      return(final_res)
    }
    state$total_reads <- state$total_reads + n_reads_in_chunk
    list(state = state, sequences=chunk, should_terminate=FALSE)
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
    read_name_table <- data.frame(read=read_names)
    payload_table <- purrr::map(filtered_res$demultiplex_res$payload,
                                as.character()) %>% list_cbind()
    
    chunk_table <- cbind(read_name_table, payload_table, barcode_table)
    if (!state$output_table_initialized) {
      append <- FALSE
      state$output_table_initialized <- TRUE
    } else {
      append <- TRUE
    }
    data.table::fwrite(x = chunk_table, file = output_table_file, append = append,
           row.names = FALSE, col.names = !append, sep = "\t", eol = "\n")
    within(state,
           {
            demultiplexed_reads <- demultiplexed_reads + nrow(barcode_matrix)
         message(glue("Processed {total_reads} reads, successfully demultiplexed {demultiplexed_reads} so far..."))
           }
    )
    state
  }
  res
}
