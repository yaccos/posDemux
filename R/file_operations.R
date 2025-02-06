#' Demultiplex combinatorial barcodes
#'
#' @param input_file The file path to the fastq sequences to be demultiplexed,
#' @param output_table_file The file path to the table
#' showing the barcode assignments of each sequences
#' @param barcode_files Named character vector, 
#' the file paths to the barcodes to be demultiplexed, arranged in the same
#' order as in the sequence annotation
#' @param segments Character vector showing the segments of the
#' sequences from 5' end to 3' end. The code applied is as follows:
#'   \itemize{
#'     \item \code{'A'}: Adapter, it is trimmed and ignored
#'     \item \code{'B'}: Barcode, used for demultiplexing
#'     \item \code{'P'}: Payload, sequence to be kept after trimming
#'     and demultiplexering (e.g. cDNA).
#'   }
#' @param segment_lengths Integer vector with the same length
#'  as \code{segments}, lengths of the segments.
#'  Up to one of the non-barcode segments can have its length
#'  set to \code{NA} which means
#'  it is considered a variadic length segment
#' @param allowed_mismatches Maximum Hamming distance from sequences
#' to barcodes, either a single integer or an integer vector with one number
#' per barcode file
#' @param trimmed_output_file The file path to the fastq files
#' containing the  sequences of the payload after adapters
#' and barcodes have been trimmed. If \code{NULL}, no such file will be created
#' @importFrom Biostrings readQualityScaledDNAStringSet
#' @importFrom Biostrings readDNAStringSet
#' @importFrom Biostrings subseq
#' @importFrom Biostrings width
#' 
#'
#' @return The output from \code{\link{filter_demultiplex_res}}. In addition,
#' the function will write to the file paths specified in \code{output_table_file}
#' and (optionally) \code{trimmed_output_file}.
#' 
#' @details
#' The barcode table will be tab separated with the first row showing the
#' barcode name.
#' names.
#' @export
#' 
#' 
file_combinatorial_demultiplex <- function(input_file,
                                      output_table_file, barcode_files,
                                      sequence_annotation, segment_lengths,
                                      allowed_mismatches=0L, 
                                      trimmed_output_file=NULL) {
  sequences <- readQualityScaledDNAStringSet(input_file,
                                             quality.scoring = "PhredQuality")
  barcodes <- map(barcode_files, readDNAStringSet, format = "fasta")
  demultiplex_res <- combinatorial_demultiplex(sequences, barcodes,
                                               sequence_annotation, segment_lengths)
  filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches)
  if (!is.null(trimmed_output_file)) {
    writeXStringSet(filtered_res$demultiplex_res$payload, trimmed_output_file)
  }
  
  write.table(filtered_res$demultiplex_res$assigned_barcodes,
              output_table_file, sep="\t", quote=FALSE, row.names=TRUE,
              col.names=TRUE)
  return(filtered_res)
}
#' @title Run pipeline for demultiplexing sequences
#' @param reverse_file The file path to the reverse reads
#'  corresponding to \code{input_file}. The contents of the reverse reads
#'  will be filtered according to the same barcode assignments as the forward reads.
#' @param reverse_output_table_file If \code{reverse_read} the file path to the table
demultiplex_pipeline <- function(input_file, reverse_file,
                                 output_table_file, barcode_files,
                                 sequence_annotation, segment_lengths,
                                 allowed_mismatches, 
                                 trimmed_output_file, reverse_output_file){
  forward_res <- file_combinatorial_demultiplex(input_file,
                                                output_table_file, barcode_files,
                                                sequence_annotation, segment_lengths,
                                                allowed_mismatches, trimmed_output_file)
  
}


