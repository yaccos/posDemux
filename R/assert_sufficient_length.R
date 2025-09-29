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
  head_msg <- glue("{n_too_short} sequence(s) are too short given the provided lengths of the segments")
  if(rlang::is_named(sequences) && n_too_short <= 100L) {
    # We print their names if they exist
    # We only print the identifiers if there 100 of them or less.
    # Printing hundreds or thousands of identifiers is not informative and 
    # in extreme cases may make R to suffer from a stack overflow
    tail_msg <- glue("These are:")
    names_too_short <- names(sequences)[!sequence_long_enough]
    collased_names <- paste(names_too_short, collapse = "\n")
    stop(glue("{head_msg}\n{tail_msg}\n{collased_names}"))
  } else {
    # Otherwise just print how many there are
    # We could in theory print their positions, but this quickly
    # gets inconvenient when performing FASTQ streaming
    stop(head_msg)
  }
}

# Internal warning function for the cases where sequences are too short
# given the segments
warn_sufficient_length <- function(sequences, minimum_length) {
  sequence_long_enough <- width(sequences) >= minimum_length
  if (all(sequence_long_enough)) {
    # All is fine, return to caller
    return(sequences)
  }
  # One or more sequences are not long enough
  n_too_short <- sum(!sequence_long_enough)
  head_msg <- glue("{n_too_short} sequence(s) in the chunk are too short given the minimum width.
                   They will be removed and ignored in the statistics.
                   ")
  if(rlang::is_named(sequences) && n_too_short <= 100L) {
    # We print their names if they exist
    # We only print the identifiers if there 100 of them or less.
    # Printing hundreds or thousands of identifiers is not informative and 
    # in extreme cases may make R to suffer from a stack overflow
    tail_msg <- glue("Their identifiers are:")
    names_too_short <- names(sequences)[!sequence_long_enough]
    collased_names <- paste(names_too_short, collapse = "\n")
    warning(glue("{head_msg}\n{tail_msg}\n{collased_names}"))
  } else {
    # Otherwise just print how many there are
    # We could in theory print their positions, but this quickly
    # gets inconvenient when performing FASTQ streaming
    warning(head_msg)
  }
  sequences[sequence_long_enough]
}
