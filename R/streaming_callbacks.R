#' @importFrom purrr list_cbind
#' @export
workflow_funs <- function(input_file, output_table_file,
                          chunk_size=1e6) {
  res <- list()
  res$state_init <- list(total_reads=0L, demultiplexed_reads=0L,
                     output_table_initialized=FALSE)
  
  res$loader <- function(state) {
    if (!state$output_table_initialized) {
      message("Initializing FASTQ stream and output table")
      state$istream <- ShortRead::FastqStreamer(input_file,
                                                        n = chunk_size)
      empty_chunk <- character() %>%  DNAStringSet()
      # We make a pass though the demultiplexer with an empty chunk in order to
      # start writing to the output barcode table
      initial_res <- list(state=state, sequences=empty_chunk,
                          should_terminate=FALSE)
      return(initial_res)
    }
    chunk  <- yield(state$istream) %>% ShortRead::sread()
    n_reads_in_chunk <- length(chunk)
    if (n_reads_in_chunk == 0L || state$output_table_initialized) {
      # The case when the initial chunk is empty is given special treatment since
      # we want to create the table regardless
      # No more reads to process, make the outer framework return
      #Clean up resources
      close(state$istream)
      state$istream <- NULL
      final_res <- list(state = state, sequences=NULL, should_terminate=TRUE)
      message(glue("Done demultiplexing {state$total_reads}"))
      return(final_res)
    }
    state$total_reads <- state$total_reads + n_reads_in_chunk
    list(state = state, sequences=NULL, should_terminate=FALSE)
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
