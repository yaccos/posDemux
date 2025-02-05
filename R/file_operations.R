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
#' @return \code{NULL}, all output will be in the form of files
#' @export
#' 
#' 
file_combinatorial_demultiplex <- function(input_file,
                                      output_table_file, barcode_files,
                                      sequence_annotation, segment_lengths,
                                      allowed_mismatches=0L,
                                      trimmed_output_file=NULL) {
  sequences <- readQualityScaledDNAStringSet(input_file, quality.scoring = "PhredQuality")
  barcodes <- map(barcode_files, readDNAStringSet, format = "fasta")
}