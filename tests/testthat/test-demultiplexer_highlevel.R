library(magrittr)
library(purrr)
library(dplyr)
library(tibble)
library(Biostrings)
library(glue)
library(ggplot2)
library(posDemux)
set.seed(29505L)

barcode_frame <- create_barcode_frame(
  segment_map = segment_map,
  segment_length = segment_length,
  n_allowed_mismatches = barcode_n_allowed_mismatches
)


possible_barcode_combinations <- barcode_frame$barcode_reference %>%
  map(names) %>% do.call(function(...)
    expand.grid(..., stringsAsFactors = FALSE), .)

prelim_freq_table <- create_preliminary_freq_table(
  possible_barcode_combinations,
  n_unique_barcodes,
  mean_reads_per_artifact,
  mean_reads_per_cell,
  rbinom_size
)

n_reads <- sum(prelim_freq_table$frequency)

# Shows the order in which the reads appear when shuffled as
# we want them to appear in random order
read_order <- sample.int(n_reads)

# The barcode combination membership of each read
# corresponding to the rows the the frequency table,
# assuming no reads are discarded
combination_membership <- rep(seq_len(n_unique_barcodes), times = prelim_freq_table$frequency) %>%
  magrittr::extract(read_order)


preliminary_assigned_barcodes <- prelim_freq_table[combination_membership, barcode_frame$name, drop =
                                                               FALSE] %>%
  as.matrix() %>%
  set_rownames(NULL)

# Ensures to shuffle barcodes
# magrittr::extract(nrow(.) %>% seq_len() %>% sample(),)

# Barcodes before mutation
preliminary_bc_stringset <- map2(barcode_frame$name, barcode_frame$barcode_reference, function(name, reference) {
  barcode <- preliminary_assigned_barcodes[, name, drop = TRUE]
  reference[barcode]
}) %>%
  set_names(barcode_frame$name)

mutation_count <- get_mutation_count(mutation_p, n_reads)

mutated_bc_stringset <- map2(mutation_count,
                             preliminary_bc_stringset,
                             mutate_barcode_stringset)

mismatches_above_threshold <- map2(mutated_bc_stringset, barcode_frame$n_allowed_mismatches, function(stringset, allowed_mismatches) {
  stringset$mismatches > allowed_mismatches
})

removed_reads <- Reduce(`|`, mismatches_above_threshold)

expected_n_removed <- removed_reads %>% sum()

expected_filtered_frequency_table <- create_filtered_freq_table(prelim_freq_table,
                                                                     combination_membership,
                                                                     removed_reads)

test_plot_generation(expected_filtered_frequency_table)

adapter_frame <- create_adapter_frame(segment_map, segment_length)
payload_frame <- create_payload_frame(segment_map,
                                      segment_length,
                                      n_reads,
                                      mean_reads_per_cell,
                                      rbinom_size)
sequences <- combine_segments(barcode_frame,
                              mutated_bc_stringset,
                              adapter_frame,
                              payload_frame)
sequence_names <- names(sequences)

expected_payload <- payload_frame$sequence %>% map(. %>%
                                                     DNAStringSet() %>%
                                                     `names<-`(sequence_names))

names(expected_payload) <- which(segment_map == "P") %>% names()

names(removed_reads) <- sequence_names

rownames(preliminary_assigned_barcodes) <- sequence_names

expected_filtered_payload <- expected_payload %>%
  map(. %>% .[!removed_reads])

demultiplex_res <- combinatorial_demultiplex(
  sequences = sequences,
  barcodes = barcode_frame$barcode_reference,
  segments = segment_map,
  segment_lengths = segment_length
)

expected_mismatches <- map(mutated_bc_stringset, "mismatches") %>%
  {
    do.call(cbind, .)
  } %>%
  set_rownames(NULL) %>%
  # We need to unclass() the rownames because they are of type glue due to the
  # tests being overly sensitive
  set_rownames(sequence_names %>% unclass())


test_that("Unfiltered payload extraction is correct", {
  compare_sequence_list(expected_payload, demultiplex_res$payload)
})

test_that("All barcodes mismatches are recorded correctly", {
  testthat::expect_true(all.equal(demultiplex_res$mismatches, expected_mismatches))
})

demultiplex_filter <- filter_demultiplex_res(demultiplex_res = demultiplex_res,
                                             allowed_mismatches = barcode_frame$n_allowed_mismatches)

expected_assigned_barcodes <- preliminary_assigned_barcodes[!removed_reads, ]

test_that("The correct reads are filtered", {
  expect_true(all.equal(demultiplex_filter$retained, !removed_reads))
})

test_that("Filtered payload extraction is correct", {
  compare_sequence_list(expected_filtered_payload,
                        demultiplex_filter$demultiplex_res$payload)
})


test_that("Barcode assignments are correct", {
  testthat::expect_true(
    all.equal(
      demultiplex_filter$demultiplex_res$assigned_barcodes,
      expected_assigned_barcodes
    )
  )
})

expected_summary_res <- create_expected_summary_res(
  possible_barcode_combinations,
  expected_filtered_frequency_table,
  n_reads,
  removed_reads,
  barcode_frame,
  mismatches_above_threshold,
  mutation_count
)

test_that("Filtering summary is correctly generated", {
  expect_equal(demultiplex_filter$summary_res, expected_summary_res)
})

freq_table <- create_freq_table(demultiplex_filter$demultiplex_res$assigned_barcodes)


test_that("Generated frequency table is correct", {
  test_freq_table(freq_table,
                       expected_filtered_frequency_table,
                       barcode_frame$name)
})
