library(magrittr)
library(purrr)
library(dplyr)
library(tibble)
library(Biostrings)
library(glue)
library(ggplot2)
library(stringi)
library(posDemux)
set.seed(29505L)
# Three barcode sets
segment_map <- c(a_1="A",bc_1="B",a_2="A",p_1="P",
                 bc_2="B",p_2="P",a_4="A",bc_3="B",p_3="P")
n_segments <- length(segment_map)
n_unique_barcodes <- 1000L
mean_reads_per_cell <- 100L
mean_reads_per_artifact <- 5L 
n_barcode_sets <- (segment_map == "B") %>% sum()
segment_length <- c(20L, 9L, 13L, 9L, 12L, NA_integer_, 18L, 5L, 14L)
barcode_frame <- tibble(segment_idx =
                          segment_map %>%
                          magrittr::equals("B") %>%
                          which()
) %>%
  mutate(length = segment_length %>% magrittr::extract(segment_idx), 
         n_allowed_mismatches = c(2L, 3L, 1L),
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
                                    ) %>% DNAStringSet()
                                  }
  ) %>% set_names(name)
  )

possible_barcode_combinations <- barcode_frame$barcode_reference %>%
  map(names) %>% do.call(function(...) expand.grid(...,stringsAsFactors = FALSE), .)
realized_barcode_combinations <- possible_barcode_combinations %>%
  slice_sample(n=n_unique_barcodes)

# We want half to the barcode combinations to simulate artifacts,
# whereas the others look more sensible to originate from cells
n_artifact_barcode_combinations <- n_unique_barcodes %/% 2L
artifact_idxs <- sample(realized_barcode_combinations %>% nrow() %>% seq_len(),
                        n_artifact_barcode_combinations)
is_artifact <- rep(FALSE, n_unique_barcodes) %>% inset(artifact_idxs, TRUE)
rbinom_size <- 10L

# How the frequency table would have
# looked if there were no errors in the barcodes
preliminary_frequency_table <-
  realized_barcode_combinations %>%  mutate(
    frequency=ifelse(is_artifact,
                     rnbinom(n_unique_barcodes,mu=mean_reads_per_artifact, size=rbinom_size),
                     rnbinom(n_unique_barcodes,mu=mean_reads_per_cell, size=rbinom_size)
    ) %>% ifelse(. == 0L, 1L, .)
  ) %>% 
  arrange(desc(frequency)) %>% 
  mutate(cumulative_frequency = cumsum(frequency), 
         fraction = frequency / sum(frequency)) %>% 
  mutate(cumulative_fraction=cumsum(fraction))

n_reads <- sum(preliminary_frequency_table$frequency)

# Shows the order in which the reads appear when shuffled as
# we want them to appear in random order
read_order <- sample.int(n_reads)


# The barcode combination membership of each read
# corresponding to the rows the the frequency table,
# assuming no reads are discarded
combination_membership <- rep(seq_len(n_unique_barcodes),
                              times = preliminary_frequency_table$frequency) %>% 
  magrittr::extract(read_order)


preliminary_assigned_barcodes <- preliminary_frequency_table[combination_membership,
                                                             barcode_frame$name, drop=FALSE] %>% 
  as.matrix() %>%
  set_rownames(NULL)

# Ensures to shuffle barcodes
# magrittr::extract(nrow(.) %>% seq_len() %>% sample(),)

# Barcodes before mutation
preliminary_bc_stringset <- map2(barcode_frame$name, barcode_frame$barcode_reference,
                                 function(name, reference) {
                                   barcode <- preliminary_assigned_barcodes[,name,drop=TRUE]
                                   reference[barcode]
                                 }
) %>%
  set_names(barcode_frame$name)

mutation_p <- list(
  # Read as: p=0.1 for 1 mismatch, p=0.04 for 2 mismatches,
  # p=0.008 for 3 mismatches (above threshold)
  bc1=c(0.1,0.04,0.008),
  bc2=c(0.2,0.05,0.02,0.005),
  bc3=c(0.05,0.01)
)



mutation_count <- mutation_p %>% map(function(p_vec) {
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

mutated_bc_stringset <- map2(mutation_count, preliminary_bc_stringset,
                             function(mutation_count, stringset) 
                             {
                               # We don't want to bother with zero mutations
                               mutation_count <- mutation_count %>% magrittr::extract(-1)
                               posDemux:::mutate_barcodes(barcodes = stringset,
                                                          mismatches = seq_along(mutation_count),
                                                          times = mutation_count)
                               
                             }
)

mismatches_above_threshold <- map2(mutated_bc_stringset,
                                   barcode_frame$n_allowed_mismatches,
                                   function(stringset, allowed_mismatches) {
                                     stringset$mismatches > allowed_mismatches 
                                   }
)

removed_reads <- Reduce(`|`, mismatches_above_threshold)

expected_n_removed <- removed_reads %>% sum()

reads_removed_per_barcode <- table(combination_membership[removed_reads])


subtraction_idxs <- names(reads_removed_per_barcode) %>% as.integer()

expected_filtered_frequency_table <- preliminary_frequency_table %>%
  mutate(frequency = frequency %>% inset(subtraction_idxs,
                                         .[subtraction_idxs] - reads_removed_per_barcode)
  ) %>% 
  filter(frequency > 0) %>% 
  arrange(desc(frequency)) %>% 
  mutate(cumulative_frequency = cumsum(frequency), 
         fraction = frequency / sum(frequency)) %>% 
  mutate(cumulative_fraction=cumsum(fraction))

test_that("Frequency and Knee plots are made without raising errors", {
  frequency_plot_file <- tempfile("frequency_plot", fileext = ".pdf")
  knee_plot_file <- tempfile("knee_plot.pdf", fileext = ".pdf")
  suppressMessages(
    # Otherwise generated annoying saving messages
    {
      expect_no_error(
        frequency_plot(expected_filtered_frequency_table) %>% ggsave(filename = frequency_plot_file, 
                                                                     plot = .)
      )
      expect_no_error(
        knee_plot(expected_filtered_frequency_table) %>% ggsave(filename = knee_plot_file, 
                                                                plot = .)
      )
    }
  )
}
)


segment_seq <- vector(mode = "list", length = n_segments)

segment_seq[barcode_frame$segment_idx] <- mutated_bc_stringset %>% map("barcodes")

sample_DNA <- function(length) {
  sample(DNA_BASES, length, replace = TRUE)
}

adapter_frame <- tibble(segment_idx = (segment_map == "A") %>% which()) %>% 
  mutate(length = segment_length[segment_idx]) %>% 
  mutate(sequence = map_chr(length, . %>% sample_DNA %>% paste0(collapse = ""))
  )


segment_seq[adapter_frame$segment_idx] <- adapter_frame$sequence

payload_frame <- tibble(segment_idx = (segment_map == "P") %>% which()) %>% 
  mutate(length = segment_length[segment_idx]) %>% 
  mutate(sequence = map(length, function(length) {
    if (is.na(length))
    {
      # Variadic segment
      # These parameters are reused outside
      # their original context because we want a similar distribution of the length
      # of the variadic segment
      length <- rnbinom(n_reads, mu=mean_reads_per_cell, size=rbinom_size) %>%
        {ifelse(. == 0L, 1L, .)}
      bulk_sampled_DNA <- sum(length) %>% sample_DNA()
      # A performance trick because using an integer vector in the subsequent
      # call to split() would result in R trying to check for unique values
      # which incur a considerable overhead in our case
      length_factor <- length %>% seq_along() %>% factor()
      res <- split(bulk_sampled_DNA, rep(length_factor, times=length)) %>% 
        map_chr(. %>% paste0(collapse = ""))
      return(res)
    }
    bulk_sampled_DNA <- sample_DNA(n_reads * length)
    split(bulk_sampled_DNA, rep(seq_len(n_reads), length)) %>%
      map_chr(.  %>% paste0(collapse = ""))
  }
  )
  )


segment_seq[payload_frame$segment_idx] <- payload_frame$sequence


sequences <- do.call(xscat, segment_seq)

sequence_names <- glue("seq_{seq_along(sequences)}")

expected_payload <- payload_frame$sequence %>% map(. %>%
                                                     DNAStringSet() %>%
                                                     `names<-`(sequence_names))
names(expected_payload) <- which(segment_map == "P") %>% names()

names(removed_reads) <- sequence_names

names(sequences) <- sequence_names

rownames(preliminary_assigned_barcodes) <- sequence_names

expected_filtered_payload <- expected_payload %>%
  map(. %>% .[!removed_reads])

demultiplex_res <- combinatorial_demultiplex(sequences = sequences,
                                             barcodes = barcode_frame$barcode_reference,
                                             segments = segment_map,
                                             segment_lengths = segment_length
)

expected_mismatches <- map(mutated_bc_stringset, "mismatches") %>%
  {do.call(cbind,.)} %>% 
  set_rownames(NULL) %>% 
  # We need to unclass() the rownames because they are of type glue due to the
  # tests being overly sensitive
  set_rownames(sequence_names %>% unclass())

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

test_that("Unfiltered payload extraction is correct",
          {
            compare_sequence_list(expected_payload, demultiplex_res$payload)
          }
)


test_that("All barcodes mismatches are recorded correctly",
          {
            testthat::expect_true(all.equal(demultiplex_res$mismatches,
                                            expected_mismatches)
            )
          }
)




demultiplex_filter <- filter_demultiplex_res(demultiplex_res = demultiplex_res,
                                             allowed_mismatches = barcode_frame$n_allowed_mismatches)



expected_assigned_barcodes <- preliminary_assigned_barcodes[!removed_reads,]

test_that("The correct reads are filtered",
          {
            expect_true(all.equal(demultiplex_filter$retained,
                                  !removed_reads))
          }
)

test_that("Filtered payload extraction is correct",
          {
            compare_sequence_list(expected_filtered_payload,
                                  demultiplex_filter$demultiplex_res$payload)
          }
)


test_that("Barcode assignments are correct",
          {
            testthat::expect_true(
              all.equal(demultiplex_filter$demultiplex_res$assigned_barcodes,
                        expected_assigned_barcodes)
            )
          }
)

expected_summary_res <- local(
  {
    n_barcode_combinations <- possible_barcode_combinations %>% nrow()
    # Note that some barcodes combination may disappear due to filtering
    # Hence, this must be taken into consideration by providing a local version
    # of n_unique_barcodes
    n_unique_barcodes <- nrow(expected_filtered_frequency_table)
    collision_lambda <- n_unique_barcodes / n_barcode_combinations
    list(
      n_reads = n_reads, 
      n_removed = removed_reads %>% sum(),
      n_barcode_sets = nrow(barcode_frame),
      n_barcode_combinations = n_barcode_combinations,
      n_unique_barcodes = n_unique_barcodes,
      collision_lambda = collision_lambda,
      expected_collisions = collision_lambda * n_unique_barcodes,
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
                                                               extract2(bc_name) %>%
                                                               magrittr::extract(seq_len(n_allowed_mismatches + 1L))
                                 )
                                 
                               )
                               
                             }
                             
      ) %>% 
        set_names(barcode_frame$name)
    ) %>% set_class("demultiplex_filter_summary")
  }
)

test_that("Filtering summary is correctly generated",
          {
            expect_equal(demultiplex_filter$summary_res, expected_summary_res)
          }
)

frequency_table <- create_frequency_table(demultiplex_filter$demultiplex_res$assigned_barcodes)

test_that("Generated frequency table is correct",
          {
            walk(c("frequency","cumulative_frequency",
                   "fraction","cumulative_fraction"), function(column) {
                     expect_equal(frequency_table[[column]],
                                  expected_filtered_frequency_table[[column]])
                   })
            # For testing whether the barcodes in the
            # frequency table are correct, we cannot just compare the columns
            # as barcode combinations with tied frequency values can appear in
            # any order. Hence, we test the barcode columns together with their
            # respective frequencies using set operations
            columns_to_test <- c(barcode_frame$name, "frequency")
            dplyr::setequal(frequency_table[columns_to_test],
                            expected_filtered_frequency_table[columns_to_test])
            
          }
)

