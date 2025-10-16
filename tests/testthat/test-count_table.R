library(magrittr)
test_that("Initialization of count table works", {
  table_handle <- posDemux:::create_count_table()
  empty_pair <- posDemux:::get_count_table(table_handle)
  expected_res <- list(encoding= integer(), frequency = integer())
  expect_equal(empty_pair, expected_res)
})

set.seed(5930)
n_chunks <- 10L

encodings <- c(52L, 1052L, 97L, 201L, 2L, 13L)
frequencies <- c(10L, 15L, 17L, 281L, 4L, 820L)
observed_encodings <- rep(encodings, times=frequencies) %>% sample()
chunk_assignments <- sample(rep(seq_len(n_chunks), length.out=length(observed_encodings)),
                            size = length(observed_encodings),
                            replace = TRUE)

growth_chunks <- split(observed_encodings, chunk_assignments)

test_that("Population of entries into count table works as expected", {
  table_handle <- posDemux:::create_count_table()
  for (chunk in growth_chunks) {
    posDemux:::add_table_entries(table_handle, chunk)
  }
  materialized_table <- posDemux:::get_count_table(table_handle)
  expect_equal(length(materialized_table$encoding), length(encodings))
  expect_equal(materialized_table$encoding %>% sort(), encodings %>% sort())
  table_order <- match(materialized_table$encoding, encodings)
  expect_equal(materialized_table$frequency, frequencies[table_order])
})
