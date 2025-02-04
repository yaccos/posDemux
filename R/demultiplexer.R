library(DNABarcodes)
library(Biostrings)
library(purrr)
library(magrittr)


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
combinatorial_demultiplex <- function(input_file,
                                      output_table_file, barcode_files,
                                      sequence_annotation, segment_lengths,
                                      allowed_mismatches=0L,
                                      trimmed_output_file=NULL) {
  sequences <- readQualityScaledDNAStringSet(input_file, quality.scoring = "PhredQuality")
  barcodes <- map(barcode_files, readDNAStringSet, format = "fasta")
  n_barcode_segments <- sum(segments == "B")
  
  if (n_barcode_segments != length(barcodes)) {
    stop("The number of barcode files does not match the number of barcode elements")
  }
  
  if(length(allowed_mismatches) != 1L && length(allowed_mismatches) != length(barcodes)) {
    stop("allowed_mismatches must be a single integer or a vector with one number per barcode file")
  }
  
  if (length(segments) != length(segment_lengths)) {
    stop("The length of segments does not match the length of segment_lengths")
  }
  
  if (is.null(names(barcode_files))) {
    stop("barcode_files must be named")
  }
  
  element_NA_idx <- which(is.na(segment_lengths))
  if (length(element_NA_idx) > 1L) {
    stop("Only one NA is allowed in segment_lengths")
  }
  
  
  
  if (length(element_NA_idx) == 1L) {
    varidic_segment_type <- segments[element_NA_idx]
    if(varidic_segment_type == "B") {
      stop("A barcode segment cannot have variadic length")
    }
    five_prime_segments <- segments[seq_len(element_NA_idx - 1L)]
    three_prime_segments <- segments[element_NA_idx +
                                                  seq_len(length(segments) - element_NA_idx)]
    five_prime_lengths <- segment_lengths[seq_len(element_NA_idx - 1L)]
    three_prime_lengths <- segment_lengths[element_NA_idx +
                                              seq_len(length(segment_lengths) - element_NA_idx)]
    n_five_prime_barcodes <- sum(five_prime_segments == "B")
    n_three_prime_barcodes <- sum(three_prime_segments == "B")
    five_prime_barcodes <- barcodes[seq_len(n_five_prime_barcodes)]
    three_prime_barcodes <- barcodes[n_five_prime_barcodes +
                                      seq_len(n_three_prime_barcodes)]
    five_prime_width <- sum(five_prime_lengths)
    three_prime_width <- sum(three_prime_lengths)
    five_prime_sequences <- subseq(sequences, start = 1L, width = five_prime_width)
    three_prime_sequences <- subseq(sequences, end = width(sequences),
                                    width = three_prime_width)
    if (length(allowed_mismatches) != 1L) {
      five_prime_allowed_mismatches <- allowed_mismatches[seq_len(n_five_prime_barcodes)]
      three_prime_allowed_mismatches <- allowed_mismatches[n_five_prime_barcodes +
                                                            seq_len(n_three_prime_barcodes)]
    } else {
      five_prime_allowed_mismatches <- rep(allowed_mismatches, n_five_prime_barcodes)
      three_prime_allowed_mismatches <- rep(allowed_mismatches, n_three_prime_barcodes)
    }
    five_prime_results <- extract_and_demultiplex(five_prime_sequences, five_prime_barcodes,
                                                  five_prime_segments, five_prime_lengths,
                                                  five_prime_allowed_mismatches)
    three_prime_results <- extract_and_demultiplex(three_prime_sequences, three_prime_barcodes,
                                                  three_prime_segments, three_prime_lengths,
                                                  three_prime_allowed_mismatches)
    
    variadic_sequence <- subseq(sequences, start = five_prime_width + 1L,
                                      end = width(sequences) - three_prime_width - 1L)
  } else {
    # No variable length payload segment
    if (length(allowed_mismatches) != 1L) {
      allowed_mismatches <- rep(allowed_mismatches, length(barcodes))
    }
    results <- extract_and_demultiplex(sequences, barcodes, segments, segment_lengths,
                                      allowed_mismatches)
  }
  
  
}

# This function assumes that the exact lengths of all segments are known
# Arguments sequences and barcodes are DNAStringSet objects
# Allowed mismatches is always an integer vector with the same length as barcodes
extract_and_demultiplex <- function(sequences, barcodes,
                                    segments, segment_lengths,
                                    allowed_mismatches) {
  barcode_widths <- imap_int(barcodes, function(barcode, name) {
    widths <- width(barcode)
    if (length(unique(widths)) > 1L) {
      stop(paste("Barcode set", name, "has variable length"))
    }
    widths[1L]
  })
  n_barcodes <- length(barcodes)
  n_segments <- length(segments)
  segment_ends <- cumsum(segment_lengths)
  barcode_segment_idxs <- which(segments == "B")
  if (!all.equal(barcode_widths, segment_lengths[barcode_segment_idxs])) {
    stop("Barcodes lengths do not match their provided segments lengths")
  }
  payload_segment_idxs <- which(segments == "P")
  barcode_segments_sequences <- map2(barcode_widths, barcode_segment_idxs,
                                    function(width, idx) {
                                      subseq(sequences,
                                            end = segment_ends[idx],
                                            width = width)
                                    })
  payload_widths <- segment_lengths[payload_segment_idxs]
  payload_sequences <- map2(payload_widths, payload_segment_idxs,
                            function(width, idx) {
                              subseq(sequences,
                                    end = segment_ends[idx],
                                    width = width)
                            })
  barcode_results <- map2(barcode_segments_sequences, barcodes, 
                          function(segment, barcode) {
                            hamming_match(segment, barcode, allowed_mismatches)
                          }
  )
  payload <- do.call(xscat, payload_sequences)
  
}

hamming_match <- function(segment, barcode, allowed_mismatches) {
  width <- width(barcode)[0L]
  distance_res <- hamming_match(segment,names(segment), barcode,
  names(barcode), width)
  
    
  
}






