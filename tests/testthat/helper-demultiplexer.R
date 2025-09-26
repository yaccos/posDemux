library(magrittr)

create_barcode_frame <- function(segment_map, segment_length, n_allowed_mismatches) {
  tibble(segment_idx =
           segment_map %>%
           magrittr::equals("B") %>%
           which()
  ) %>%
    mutate(length = segment_length %>% magrittr::extract(segment_idx), 
           n_allowed_mismatches = n_allowed_mismatches,
           name = glue("bc{1L:3L}")) %>% 
    mutate(barcode_reference = pmap(list(name, length, n_allowed_mismatches),
                                    function(name, length, n_allowed_mismatches)
                                    {
                                      DNABarcodes::create.dnabarcodes(n=length,
                                                                      # In order to uniquely error-check a barcode
                                                                      # the distance between the barcodes have to be
                                                                      # more than double the redundancy
                                                                      # In order to also avoid misclassifications
                                                                      # of barcodes with one more mismatch than the redundancy,
                                                                      # we have to increase the distance by one more nucleotide
                                                                      dist = n_allowed_mismatches * 2L + 2L,
                                                                      metric = "hamming",
                                                                      heuristic = "conway",
                                                                      # While this kind of filtering might be useful for real
                                                                      # DNA barcodes, it serves no purpose when only the computer
                                                                      # has to deal with it
                                                                      filter.triplets = FALSE,
                                                                      filter.gc = FALSE,
                                                                      filter.self_complementary = FALSE
                                      ) %>% set_names(
                                        {glue("{name}_{seq_len(length(.))}")}
                                      ) %>% Biostrings::DNAStringSet()
                                    }
    ) %>% set_names(name)
    )
}

create_preliminary_frequency_table <- function(possible_barcode_combinations,
                                               n_unique_barcodes,
                                               mean_reads_per_artifact,
                                               mean_reads_per_cell,
                                               rbinom_size) {
  realized_barcode_combinations <- possible_barcode_combinations %>%
    slice_sample(n = n_unique_barcodes)
  
  # We want half to the barcode combinations to simulate artifacts,
  # whereas the others look more sensible to originate from cells
  n_artifact_barcode_combinations <- n_unique_barcodes %/% 2L
  artifact_idxs <- sample(
    realized_barcode_combinations %>% nrow() %>% seq_len(),
    n_artifact_barcode_combinations
  )
  is_artifact <- rep(FALSE, n_unique_barcodes) %>% magrittr::inset(artifact_idxs, TRUE)
  
  
  # How the frequency table would have
  # looked if there were no errors in the barcodes
  preliminary_frequency_table <-
    realized_barcode_combinations %>%  mutate(frequency = ifelse(
      is_artifact,
      rnbinom(n_unique_barcodes, mu = mean_reads_per_artifact, size =
                rbinom_size),
      rnbinom(n_unique_barcodes, mu = mean_reads_per_cell, size =
                rbinom_size)
    ) %>% ifelse(. == 0L, 1L, .)) %>%
    complete_frequency_table()
  preliminary_frequency_table
}

get_mutation_count <- function(mutation_p, n_reads) {
  mutation_p %>% map(function(p_vec) {
    # We must include the probability of zero mutations to sample correctly
    # from the multinomial distribution
    augmented_p_vec <- c(1L-sum(p_vec), p_vec)
    # Computes the number of reads having 0,1,..,length(p_vec) mutations
    # in the barcode
    # for each realized barcode combination
    rmultinom(n = 1L, size = n_reads, prob = augmented_p_vec) %>%
      as.vector()
  }
  )
}

test_plot_generation <- function(frequency_table) {
  test_that("Frequency and Knee plots are made without raising errors", {
    frequency_plot_file <- tempfile("frequency_plot", fileext = ".pdf")
    knee_plot_file <- tempfile("knee_plot.pdf", fileext = ".pdf")
    suppressMessages(
      # Otherwise generated annoying saving messages
      {
        expect_no_error(
          frequency_plot(frequency_table) %>% ggsave(filename = frequency_plot_file, 
                                                     plot = .)
        )
        expect_no_error(
          knee_plot(frequency_table) %>% ggsave(filename = knee_plot_file, 
                                                plot = .)
        )
      }
    )
  }
  )
}

compare_sequence_list <- function(actual, expect) {
  
  one_has_names <- function(x, y){
    !is.null(names(x)) || !is.null(names(y))
  }
  
  if (one_has_names(actual, expect)) {
    expect_true(all(names(actual) == names(expect)))
  }
  
  map2(actual, expect, function(actual, expect) {
    if (one_has_names(actual, expect)){
      expect_true(all(names(actual) ==  names(expect)))
    }
    expect_true(all(actual ==  expect))
  }
  )
}

create_filtered_frequency_table <- function(frequency_table,
                                            combination_membership, removed_reads) {
  reads_removed_per_barcode <- table(combination_membership[removed_reads])
  
  subtraction_idxs <- names(reads_removed_per_barcode) %>% as.integer()
  
  filtered_frequency_table <- frequency_table %>%
    mutate(frequency = frequency %>% magrittr::inset(subtraction_idxs,
                                           .[subtraction_idxs] - reads_removed_per_barcode)
    ) %>% 
    filter(frequency > 0) %>%
    complete_frequency_table()
  filtered_frequency_table
}

sample_DNA <- function(length) {
  sample(Biostrings::DNA_BASES, length, replace = TRUE)
}

complete_frequency_table <- . %>%
  arrange(desc(frequency)) %>% 
  mutate(cumulative_frequency = cumsum(frequency), 
         fraction = frequency / sum(frequency)) %>% 
  mutate(cumulative_fraction=cumsum(fraction))

create_payload_frame <- function(segment_map, segment_length, n_reads, mean_length, rbinom_size) {
  tibble(segment_idx = (segment_map == "P") %>% which()) %>% 
    mutate(length = segment_length[segment_idx]) %>% 
    mutate(sequence = map(length, function(length) {
      if (is.na(length))
      {
        # Variadic segment
        # These parameters are reused outside
        # their original context because we want a similar distribution of the length
        # of the variadic segment
        length <- rnbinom(n_reads, mu=mean_length, size=rbinom_size) %>%
          {ifelse(. == 0L, 1L, .)}
        bulk_sampled_DNA <- sum(length) %>% sample_DNA()
        # A performance trick because using an integer vector in the subsequent
        # call to split() would result in R trying to check for unique values
        # which incur a considerable overhead in our case
        length_factor <- length %>% seq_along() %>% factor()
        res <- split(bulk_sampled_DNA, rep(length_factor, times=length)) %>% 
          purrr::map_chr(. %>% paste0(collapse = ""))
        return(res)
      }
      bulk_sampled_DNA <- sample_DNA(n_reads * length)
      split(bulk_sampled_DNA, rep(seq_len(n_reads), length)) %>%
        purrr::map_chr(.  %>% paste0(collapse = ""))
    }
    )
    )
}

mutate_barcode_stringset <- function(mutation_count, stringset) 
{
  # We don't want to bother with zero mutations
  mutation_count <- mutation_count %>% magrittr::extract(-1)
  posDemux:::mutate_barcodes(barcodes = stringset,
                             mismatches = seq_along(mutation_count),
                             times = mutation_count)
}

create_adapter_frame <- function(segment_map, segment_length) {
  adapter_frame <- tibble(segment_idx = (segment_map == "A") %>% which()) %>% 
    mutate(length = segment_length[segment_idx]) %>% 
    mutate(sequence = purrr::map_chr(length, . %>% sample_DNA %>% paste0(collapse = ""))
    )
}

combine_segments <- function(barcode_frame, mutated_bc_stringset,
                             adapter_frame, payload_frame) {
  segment_seq <- vector(mode = "list", length = n_segments)
  segment_seq[barcode_frame$segment_idx] <- mutated_bc_stringset %>% map("barcodes")
  segment_seq[adapter_frame$segment_idx] <- adapter_frame$sequence
  segment_seq[payload_frame$segment_idx] <- payload_frame$sequence
  sequences <- do.call(xscat, segment_seq)
  sequence_names <- glue("seq_{seq_along(sequences)}")
  names(sequences) <- sequence_names
  sequences
}

create_expected_summary_res <- function(possible_barcode_combinations, expected_filtered_frequency_table, n_reads,
                                        removed_reads, barcode_frame, mismatches_above_threshold, mutation_count) {
  n_barcode_combinations <- possible_barcode_combinations %>% nrow()
  # Note that some barcodes combination may disappear due to filtering
  # Hence, this must be taken into consideration by providing a local version
  # of n_unique_barcodes
  n_unique_barcodes <- nrow(expected_filtered_frequency_table)
  observed_collision_lambda <- n_unique_barcodes / n_barcode_combinations
  n_estimated_features <- posDemux:::poisson_correct_n(n_barcode_combinations,
                                                    n_unique_barcodes)
  corrected_collision_lambda <- n_estimated_features / n_barcode_combinations
  expected_collisions <- 
    posDemux:::poisson_estimate_collisions(n_barcode_combinations,
                                           corrected_collision_lambda) 
  list(
    n_reads = n_reads, 
    n_removed = removed_reads %>% sum(),
    n_barcode_sets = nrow(barcode_frame),
    n_barcode_combinations = n_barcode_combinations,
    n_unique_barcodes = n_unique_barcodes,
    n_estimated_features = n_estimated_features,
    observed_collision_lambda = observed_collision_lambda,
    corrected_collision_lambda = corrected_collision_lambda,
    expected_collisions = expected_collisions,
    barcode_summary = imap(barcode_frame$name,
                           function(bc_name, i)
                           {
                             n_allowed_mismatches <- barcode_frame$n_allowed_mismatches[i]
                             names(n_allowed_mismatches) <- bc_name
                             list(
                               width = barcode_frame$length[i],
                               n_barcodes = barcode_frame$barcode_reference[[i]] %>% length(),
                               n_allowed_mismatches = n_allowed_mismatches,
                               n_removed = mismatches_above_threshold[[bc_name]] %>% sum(),
                               mismatch_frame = data.frame(n_mismatches = 0L:n_allowed_mismatches,
                                                           frequency = mutation_count %>%
                                                             magrittr::extract2(bc_name) %>%
                                                             magrittr::extract(seq_len(n_allowed_mismatches + 1L))
                               )
                               
                             )
                             
                           }
                           
    ) %>% 
      magrittr::set_names(barcode_frame$name)
  ) %>% magrittr::set_class("demultiplex_filter_summary")
}

test_frequency_table <- function(actual, expected, barcode_names) {
    purrr::walk(c("frequency","cumulative_frequency",
                  "fraction","cumulative_fraction"), function(column) {
                    expect_equal(actual[[column]],
                                 expected[[column]])
                  })
    # For testing whether the barcodes in the
    # frequency table are correct, we cannot just compare the columns
    # as barcode combinations with tied frequency values can appear in
    # any order. Hence, we test the barcode columns together with their
    # respective frequencies using set operations
    columns_to_test <- c(barcode_names, "frequency")
    expect_true(dplyr::setequal(actual[columns_to_test],
                    expected[columns_to_test]))
  }
