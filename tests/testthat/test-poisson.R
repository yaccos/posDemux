# Internal testing of functions for estimating number of cells and collisions
library(magrittr)
N <- 96^3 %>% as.integer()
ns <- c(10L, 30L, 100L, 1000L, 10000L, 50000L, 100000L, 200000L, 500000L)
lambdas <- ns / N
f <- \(n) N * (1 - exp(-n/N))
g <- \(n) N * (1 - exp(-n/N) - n / N* exp(-n/ N))
g <- \(n) N * (1 - exp(-n/N) - n / N* exp(-n/ N))
test_that("The Poisson method is close to the ratio method for low n",
          {
          local({
            n <- 5000L
            poisson_estimate <- posDemux:::poisson_estimate_collisions(N,
                                                                       n / N)
            ratio_estimate <- n^2/ N / 2
            expect_lt(abs(poisson_estimate-ratio_estimate), 1.0)
          })
          }
          )

# Estimate n_real from n_obs, using a numeric method
correct_n <- function(n_obs, tol = 0.01) {
  f <- \(n) N * (1 - exp(-n/N)) - n_obs
  f_prime <- \(n) exp(-n/N)
  n <- n_obs
  # Uses Newton's method to correct the number of cells present
  f_n <- f(n)
  while(abs(f_n) |> (\(err) err > tol)() |> any()) {
    n <- n - f_n / f_prime(n)
    f_n <- f(n)
  1}
  n
}

test_that("Correcting the number of features work", 
          {
            local({
            n_obs <- N * (1- exp(-lambdas))
            n_corrected <- posDemux:::poisson_correct_n(N, n_obs)
            expect_equal(object = n_corrected, ns)
          })
          }
)

test_that("The iterative method for finding number
          of features agrees with the algebraic one", 
          {
          local({
            tol <- 10^(-5)
            n_obs <- N * (1- exp(-lambdas))
            n_corrected <- correct_n(n_obs, tol)
            expect_true(abs(ns-n_corrected) %>%
                          magrittr::is_less_than(2*tol) %>%
                          all()
                        )
          })
          }
)
