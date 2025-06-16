test_that("Initialization of dynamic encoding vector works", {
  vector_handle <- posDemux:::create_encoding_vector()
  empty_vector <- posDemux:::get_encoding_vector(vector_handle)
  expect_equal(empty_vector, integer())
})

growth_chunks <- list(c(3L, -58L, 42L, 48L, -46L, 11L), c(-34L, 88L, 33L),
                       c(-30L), numeric(), c(93L, 56L))

test_that("Population of numbers into encoding vector works as expected", {
  vector_handle <- posDemux:::create_encoding_vector()
  for (chunk in growth_chunks) {
    posDemux:::grow_encoding_vector(vector_handle, chunk)
  }
  materialized_vector <- posDemux:::get_encoding_vector(vector_handle)
  expected_vector <- unlist(growth_chunks)
  expect_equal(materialized_vector, expected_vector)
})
