library(magrittr)
library(purrr)
library(dplyr)
library(tibble)
library(Biostrings)
library(glue)
set.seed(29505L)
# Three barcode sets
segment_map <- c("A","B","A","P","B","P","A","B","P")
# n_reads <- as.integer(1e6)
n_unique_barcodes <- 10000L
mean_reads_per_cell <- 100L
mean_reads_per_artifact <- 5L 
n_barcode_sets <- (segment_map == "B") %>% sum()
segment_length <- c(20L, 8L, 13L, 9L, 9L, NA_integer_, 18L, 4L, 14L)
barcode_frame <- tibble(segment_idx =
                          segment_map %>%
                          equals("B") %>%
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
  map(names) %>% do.call(expand.grid, .)
realized_barcode_combinations <- possible_barcode_combinations %>%
  slice_sample(n=n_unique_barcodes)

# We want half to the barcode combinations to simulate artifacts,
# whereas the others look more sensible to originate from cells
n_artifact_barcode_combinations <- as.integer(n_unique_barcodes / 2)
artifact_idxs <- sample(realized_barcode_combinations %>% nrow() %>% seq_len(),
                        n_artifact_barcode_combinations)
is_artifact <- rep(FALSE, n_unique_barcodes) %>% inset(artifact_idxs, TRUE)
# realized_barcode_combinations$frequency <- 
#   ifelse(is_artifact,
#          rpois(n_unique_barcodes,mean_reads_per_artifact),
#          rpois(n_unique_barcodes,mean_reads_per_cell)
#          )
rbinom_size <- 10L
expected_frequency_table <-
  realized_barcode_combinations %>%  mutate(
    frequency=ifelse(is_artifact,
         rnbinom(n_unique_barcodes,mu=mean_reads_per_artifact,size=rbinom_size),
         rnbinom(n_unique_barcodes,mu=mean_reads_per_cell,size=rbinom_size)
         ) %>% ifelse(. == 0L, 1L, .)
    ) %>% 
      arrange(desc(frequency)) %>% 
      mutate(cumulative_frequency = cumsum(frequency), 
             fraction = frequency / sum(frequency)) %>% 
      mutate(cumulative_fraction=cumsum(fraction))
# When used interactively: Validate plots look as expected
# posDemux::frequency_plot(expected_frequency_table)
# posDemux::knee_plot(expected_frequency_table)

