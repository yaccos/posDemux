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
expected_fraction <- sorted_expected_frequency / sum(expected_frequency)
expected_cumulative_fraction <- cumsum(expected_fraction)
expected_unique_mat <- unique_rows %>%
  extract(., expected_order) %>%
  do.call(rbind, .)

frequency_table <- create_frequency_table(mat)

test_that("Frequency table is correct for unnamed integer matrices", {
  expect_equal(nrow(frequency_table),n_expected_combinations)
  frequency <- frequency_table$frequency
  expect_equal(frequency,sorted_expected_frequency)
  cumulative_frequency <- frequency_table$cumulative_frequency
  expect_equal(cumulative_frequency, expected_cumulative_frequency)
  # We must unname the matrix because the data frame methods
  # in create_frequency_table makes up column names if there are none
  reconstructed_unique_mat <- frequency_table %>%
    select(-frequency, -cumulative_frequency, -fraction, -cumulative_fraction) %>%
     as.matrix() %>% 
    unname()
  expect_equal(reconstructed_unique_mat, expected_unique_mat)
})

test_that("Frequency table is correct for named integer matrices", {
  named_unique_rows <- unique_rows %>% 
    map(. %>% setNames(letters[seq_along(.)]))
  named_mat <- n_expected_combinations %>%
    seq_len() %>%
    rep(times=expected_frequency) %>% 
    sample() %>% 
    extract(named_unique_rows, .) %>% 
    do.call(rbind, .)
  named_frequency_table <- create_frequency_table(named_mat)
  expected_unique_mat <- named_unique_rows %>%
    extract(., expected_order) %>%
    do.call(rbind, .)
  expect_equal(nrow(named_frequency_table),n_expected_combinations)
  frequency <- named_frequency_table$frequency
  expect_equal(named_frequency_table$frequency, sorted_expected_frequency)
  cumulative_frequency <- named_frequency_table$cumulative_frequency
  expect_equal(cumulative_frequency, expected_cumulative_frequency)
  reconstructed_unique_mat <- named_frequency_table %>%
    select(-frequency, -cumulative_frequency, -fraction, -cumulative_fraction) %>%
    as.matrix()
  expect_equal(reconstructed_unique_mat, expected_unique_mat)
})
