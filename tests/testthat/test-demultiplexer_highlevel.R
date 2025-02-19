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
segment_map <- c("A","B","A","P","B","P","A","B","P")
n_segments <- length(segment_map)
n_unique_barcodes <- 10000L
mean_reads_per_cell <- 100L
mean_reads_per_artifact <- 5L 
n_barcode_sets <- (segment_map == "B") %>% sum()
segment_length <- c(20L, 8L, 13L, 9L, 9L, NA_integer_, 18L, 4L, 14L)
barcode_frame <- tibble(segment_idx =
                          segment_map %>%
                          magrittr::equals("B") %>%
                          which()
) %>%
  mutate(length = segment_length %>% extract(segment_idx), 
         n_allowed_mismatches = c(4L, 5L, 1L),
         name = glue("bc{1L:3L}")) %>% 
  mutate(barcode_reference = pmap(list(name, length, n_allowed_mismatches),
                                  function(name, length, n_allowed_mismatches)
                                  {
                                    DNABarcodes::create.dnabarcodes(n=length,
                                                                    dist = n_allowed_mismatches - 1L ,
                                                                    metric = "hamming",
                                                                    heuristic = "conway"
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
n_artifact_barcode_combinations <- as.integer(n_unique_barcodes / 2)
artifact_idxs <- sample(realized_barcode_combinations %>% nrow() %>% seq_len(),
                        n_artifact_barcode_combinations)
is_artifact <- rep(FALSE, n_unique_barcodes) %>% inset(artifact_idxs, TRUE)
rbinom_size <- 10L
expected_frequency_table <-
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

n_reads <- sum(expected_frequency_table$frequency)

test_that("Frequency and Knee plots are made without raising errors", {
  frequency_plot_file <- tempfile("frequency_plot", fileext = ".pdf")
  knee_plot_file <- tempfile("knee_plot.pdf", fileext = ".pdf")
  expect_no_error(
    frequency_plot(expected_frequency_table) %>% ggsave(filename = frequency_plot_file, 
                                                        plot = .)
  )
  expect_no_error(
    knee_plot(expected_frequency_table) %>% ggsave(filename = knee_plot_file, 
                                                   plot = .)
  )
}
)

expected_assigned_barcodes <- map2(barcode_frame$name, barcode_frame$barcode_reference,
                                   function(name, reference) {
                                     barcode <- expected_frequency_table[[name]]
                                     barcode_multiplicity <- expected_frequency_table$frequency
                                     barcode_with_multiplicity <- rep(barcode, barcode_multiplicity)
                                     reference[barcode_with_multiplicity] %>%
                                       names()
                                   }
) %>%
  {do.call(cbind, .)} %>% 
  set_colnames(barcode_frame$name) %>% 
  # Ensures to shuffle barcodes
  extract(nrow(.) %>% seq_len() %>% sample(),)

bc_stringset <- map2(barcode_frame$name, barcode_frame$barcode_reference,
                    function(name, reference) {
                      barcode <- expected_frequency_table[[name]]
                      barcode_multiplicity <- expected_frequency_table$frequency
                      barcode_with_multiplicity <- rep(barcode, barcode_multiplicity)
                      expected_assigned_barcodes %>% extract(, name) %>% 
                        {extract(reference, .)}
                    }
                    ) %>% 
  set_names(barcode_frame$name)

segment_seq <- vector(mode = "list", length = n_segments)

segment_seq[barcode_frame$segment_idx] <- bc_stringset

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
      map_chr(. %>% paste0(collapse = ""))
  }
  )
  )

  
segment_seq[payload_frame$segment_idx] <- payload_frame$sequence

expected_payload <- do.call(xscat, payload_frame$sequence) %>% DNAStringSet()

sequences <- do.call(xscat, segment_seq)

demultiplex_res <- combinatorial_demultiplex(sequences = sequences,
                                                       barcodes = barcode_frame$barcode_reference,
                                                       segments = segment_map,
                                                       segment_lengths = segment_length)

test_that("Barcode assignments are correct",
          {
          testthat::expect_true(
            all.equal(demultiplex_res$assigned_barcodes,
                                          expected_assigned_barcodes)
          )
          }
)

expected_mismatches <- array(0L, dim = dim(expected_assigned_barcodes), 
                             dimnames = dimnames(expected_assigned_barcodes))

test_that("All barcodes were assigned without any mismatches",
          {
          testthat::expect_true(all.equal(demultiplex_res$mismatches,
                                          expected_mismatches)
          )
          }
)

demultiplex_filter <- filter_demultiplex_res(demultiplex_res = demultiplex_res,
                       allowed_mismatches = barcode_frame$n_allowed_mismatches)

frequency_table <- create_frequency_table(demultiplex_filter$demultiplex_res$assigned_barcodes)




test_that("Generated frequency table is correct",
          {
            walk(c("frequency","cumulative_frequency",
                   "fraction","cumulative_fraction"), function(column) {
                     expect_equal(frequency_table[[column]],
                                  expected_frequency_table[[column]])
                   })
            # For testing whether the barcodes in the
            # frequency table are correct, we cannot just compare the columns
            # as barcode combinations with tied frequency values can appear in
            # any order. Hence, we test the barcode columns together with their
            # respective frequencies using set operations
            columns_to_test <- c(barcode_frame$name, "frequency")
            dplyr::setequal(frequency_table[columns_to_test],
                            expected_frequency_table[columns_to_test])
            
          }
          )


