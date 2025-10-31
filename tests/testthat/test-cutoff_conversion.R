# The result of the test is not really dependent on the random seed,
# as the randomness is only used for shuffling the rows of the matrix.
set.seed(5030)
library(magrittr)
library(dplyr)
library(purrr)
unique_rows <- list(c(1L, 5L, 3L), c(4L, 5L, 6L), c(1L, 2L, 3L), c(2L, 3L, 3L))
n_expected_combinations <- length(unique_rows)
expected_frequency <- c(4L,2L,3L,3L)
expected_order <- order(expected_frequency, decreasing = TRUE)
sorted_expected_frequency <- sort(expected_frequency, decreasing = TRUE)
expected_cumulative_frequency <- cumsum(sorted_expected_frequency)
mat <- n_expected_combinations %>%
  seq_len() %>%
  rep(times=expected_frequency) %>% 
  sample() %>% 
  extract(unique_rows, .) %>% 
  do.call(rbind, .)
freq_table <- create_freq_table(mat)
bc_cutoff <- seq_len(n_expected_combinations) %>% c(0L, .)
expected_frequency_cutoff <- c(5L, 4L, 3L, 3L, 2L)



test_that("Barcode to frequency conversion works",
          {
          bc_cutoff <- seq_len(n_expected_combinations) %>% c(0L, .)
          frequency_cutoff <- bc_to_freq_cutoff(freq_table,
                                                     bc_cutoff)
          expected_frequency_cutoff <- c(5L, 4L, 3L, 3L, 2L)
          expect_equal(frequency_cutoff,
                       expected_frequency_cutoff)
          }
          )

test_that("Frequency to barcode conversion works",
          {
            # This is intended to be a decreasing sequence and actually one
            # of the cases where this behavior of the colon operator is useful
            frequency_cutoff <- 5L:0L
            bc_cutoff <- freq_to_bc_cutoff(freq_table,
                                                       frequency_cutoff)
            expected_bc_cutoff <- c(0L, 1L, 3L, 4L, 4L, 4L)
            expect_equal(bc_cutoff,
                         expected_bc_cutoff)
          }
)
