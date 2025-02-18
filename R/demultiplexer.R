#' Combinatorial demultiplexer
#'
#' @param sequences A \code{\link{XStringSet}} object, the sequences to be demultiplexed
#' @param barcodes A list of \code{\link{XStringSet}} objects, the barcodes
#' to be used for demultiplexing. All of the barcodes in each \code{\link{XStringSet}} must
#' have the same length as specified by the \code{segment_lengths} argument and be named.
#' @param segments Character vector showing the segments of the
#' sequences from 5' end to 3' end. The code applied is as follows:
#'   \itemize{
#'     \item \code{'A'}: Adapter, it is trimmed and ignored
#'     \item \code{'B'}: Barcode, used for demultiplexing
#'     \item \code{'P'}: Payload, sequence to be kept after trimming
#'     and demultiplexering (e.g. cDNA or UMI).
#'   }
#' @param segment_lengths Integer vector with the same length
#'  as \code{segments}, lengths of the segments.
#'  Up to one of the non-barcode segments can have its length
#'  set to \code{NA} which means
#'  it is considered a variadic length segment
#' @description
#' This function performs combinatorial demultiplexing and trimming
#' of sequences. The sequences have a structure as defined by the argument
#' \code{segments} with corresponding lengths in \code{segment_lengths}. 
#' As trimming and extraction can take place from either end, a single 
#' middle segment can be variadic in length. There are three types of segments:
#' Adapter, Barcode and Payload. The adapter is trimmed and ignored, the barcode
#' is used for demultiplexing, and the payload is kept after trimming and demultiplexing
#' and returned from the function. If there are multiple payload segments, they
#' concatented in the output sequences. The barcodes can be positioned at 
#' either end of the sequences,
#' but no barcode can (for obvious reasons) be variadic in length.
#' 
#' 
#' 
#' @importFrom Biostrings subseq width xscat
#'
#' @return A list with the following elements:
#'  \itemize{
#'  \item \code{assigned_barcodes}: A \code{character} matrix with
#'  the names of the assigned barcodes as elements. The rows correspond to the
#'  sequences and the columns to the barcode segments
#'  \item \code{mismatches}: An \code{integer} matrix with the number of mismatches
#'  between the assigned barcodes and the sequences. The rows correspond to the
#'  sequences and the columns to the barcode segments.
#'  \item \code{payload}: A \code{\link{XStringSet}} object with the payload sequences
#'  }
#'  
#' @details
#' If there are two barcodes both having the minimum number of mismatches
#' the first one will be selected. If is therefore important to choose the
#' error tolerance to be equal or less than the redunancy of the barcodes.
#' @export
#' 
#' 
combinatorial_demultiplex <- function(sequences, barcodes,
                                      segments, segment_lengths) {
  assert_that(is(sequences, "XStringSet"),
              msg="The argument sequences must be an XStringSet object")
  n_barcode_segments <- sum(segments == "B")
  n_segments <- length(segments)
  
  # The last line rewritten with assert_that
  assert_that(n_barcode_segments == length(barcodes),
              msg="The provided number of barcode segments in argument barcodes
         does not match the number of barcode segments in argument segments")
  
  assert_that(length(segments) == length(segment_lengths),
              msg="The length of segments does not match the length of segment_lengths")
  
  
  iwalk(barcodes, function(barcode, name) {
    assert_that(is(barcode, "XStringSet"),
                msg=glue("The barcodes of segment {name} must be an XStringSet object"))
    assert_that(!is.null(names(barcode)),
                msg=glue("The barcodes of segment {name} are not named"))
    assert_that(length(unique(width(barcode))) <= 1L,
                msg=paste("The barcodes of segment {name} have variable length"))
  }
  
  )
  
  element_NA_idx <- which(is.na(segment_lengths))
  
  assert_that(length(element_NA_idx) <= 1L,
              msg="Only one NA is allowed in segment_lengths")
  
  
  
  if (length(element_NA_idx) == 1L) {
    varidic_segment_type <- segments[element_NA_idx]
    assert_that(varidic_segment_type != "B",
                msg="A barcode segment cannot have variadic length")
    n_five_prime_segments <- element_NA_idx - 1L
    n_three_prime_segments <- n_segments - 1L - n_five_prime_segments
    five_prime_segments <- segments[seq_len(n_five_prime_segments)]
    three_prime_segments <- segments[seq.int(element_NA_idx + 1L, length.out = n_three_prime_segments)]
    five_prime_lengths <- segment_lengths[seq_len(n_five_prime_segments)]
    three_prime_lengths <- segment_lengths[seq.int(element_NA_idx + 1L, length.out = n_three_prime_segments)]
    n_five_prime_barcodes <- sum(five_prime_segments == "B")
    n_three_prime_barcodes <- sum(three_prime_segments == "B")
    five_prime_barcodes <- barcodes[seq_len(n_five_prime_barcodes)]
    three_prime_barcodes <- barcodes[seq.int(n_five_prime_barcodes + 1L, length.out = n_three_prime_barcodes)]
    five_prime_width <- sum(five_prime_lengths)
    three_prime_width <- sum(three_prime_lengths)
    five_prime_sequences <- subseq(sequences, start = 1L, width = five_prime_width)
    three_prime_sequences <- subseq(sequences, end = width(sequences),
                                    width = three_prime_width)
    five_prime_results <- extract_and_demultiplex(five_prime_sequences, five_prime_barcodes,
                                                  five_prime_segments, five_prime_lengths)
    three_prime_results <- extract_and_demultiplex(three_prime_sequences, three_prime_barcodes,
                                                   three_prime_segments, three_prime_lengths)
    assigned_barcodes <- cbind(five_prime_results$assigned_barcode,
                               three_prime_results$assigned_barcode)
    mismatches <- cbind(five_prime_results$mismatches,
                        three_prime_results$mismatches)
    if (varidic_segment_type == "P") {
      variadic_sequence <- subseq(sequences, start = five_prime_width + 1L,
                                  end = width(sequences) - three_prime_width - 1L)
      payload <- xscat(five_prime_results$payload, variadic_sequence,
                       three_prime_results$payload)
    } else {
      payload <- xscat(five_prime_results$payload, three_prime_results$payload)
    }
    return(list(assigned_barcodes = assigned_barcodes,
                mismatches = mismatches,
                payload = payload))
    
    
  } else {
    # No variable length payload segment
    results <- extract_and_demultiplex(sequences, barcodes, segments, segment_lengths)
  }
  return(results)
}

# This function assumes that the exact lengths of all segments are known
# Arguments sequences and barcodes are DNAStringSet objects
extract_and_demultiplex <- function(sequences, barcodes,
                                    segments, segment_lengths) {
  barcode_widths <- imap_int(barcodes, function(barcode, name) {
    widths <- width(barcode)
    # The last line rewritten with assert_that
    assert_that(length(unique(widths)) <= 1L,
                msg=glue("Barcode set {name} has variable length"))
    widths[1L]
  })
  n_barcodes <- length(barcodes)
  n_segments <- length(segments)
  segment_ends <- cumsum(segment_lengths)
  barcode_segment_idxs <- which(segments == "B")
  assert_that(all.equal(barcode_widths, segment_lengths[barcode_segment_idxs],
                        check.attributes = FALSE) %>% 
                isTRUE(),
              msg="Barcodes lengths do not match their provided segments lengths")
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
  barcode_results <- pmap(list(barcode_segments_sequences, barcodes, barcode_widths), 
                          function(segment, barcode, width) {
                            hamming_match(segment, names(segment), barcode,
                                          names(barcode), width)
                          }
  )
  assigned_barcode <- do.call(cbind, map(barcode_results, "assigned_barcode"))
  mismatches <- do.call(cbind, map(barcode_results, "mismatches"))
  payload <- do.call(xscat, payload_sequences)
  list(assigned_barcode = assigned_barcode,
              mismatches = mismatches,
              payload = payload) %>% return()
}

