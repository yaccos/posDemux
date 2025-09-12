# Internal error-handling function for the cases where sequences are too short
# given the segments
assert_sufficient_length <- function(sequences, minimum_length) {
  sequence_long_enough <- width(sequences) >= minimum_length
  if (all(sequence_long_enough)) {
    # All is fine, return to caller
    return()
  }
  # One or more sequences are not long enough
  n_too_short <- sum(!sequence_long_enough)
  head_msg <- glue("{n_too_short} sequence(s) are too short given the provided lengths of the segments\n")
  rlang::is_named(sequences){
    # We print their names if they exist
    tail_msg <- glue("These are:")
    names_too_short <- names(sequences)
    collased_names <- paste(names_too_short, collapse = "\n")
    stop(head_msg, tail_msg, collased_names)
  } else {
    # Otherwise just print how many there are
    # We could in theory print their positions, but this quickly
    # gets inconvenient when performing FASTQ streaming
    stop(head_msg)
  }
}
