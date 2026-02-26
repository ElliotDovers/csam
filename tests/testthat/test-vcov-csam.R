# tests/testthat/test-vcov-csam.R
library(testthat)

# helper to simulate tiny data and fit csam
make_and_fit <- function(n = 10, s = 10, p = 1, g = 2, d = 1, family = poisson(), psi1 = 1, psi2 = 1) {
  set.seed(123)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  # simulate simple Poisson counts with small mean
  Y <- matrix(rpois(n * s, lambda = 3), nrow = n, ncol = s)
  # Fit csam; use small max_iter to keep tests fast
  fit <- csam(Y = Y, X = X, g = g, d = d, family = family,
              psi1 = psi1, psi2 = psi2, max_iter = 5, trace = TRUE,
              maxit_step1 = 1, maxit_step2 = 1, maxit_step3 = 1, verbose = FALSE)
  fit
}

# small tolerance for eigen negative due to numerical error
tol_eig <- 1e-6

test_that("vcov.csam works for g = 1 (no free pi)", {
  fit <- make_and_fit(g = 1)
  expect_silent({
    V_naive <- vcov(fit, method = "naive")
    V_louis <- vcov(fit, method = "louis")
    V_oakes <- vcov(fit, method = "oakes")
    V_sand  <- vcov(fit, method = "sandwich")
  })
  # basic checks
  for (V in list(V_naive, V_louis, V_oakes, V_sand)) {
    expect_true(is.matrix(V))
    expect_true(isSymmetric(V, tol = 1e-8))
    expect_true(all(is.finite(V)))
    # near positive semidefinite
    eigs <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
    expect_true(min(eigs) > -tol_eig)
  }
})

test_that("vcov.csam works for g = 2 (single free pi parameter)", {
  fit <- make_and_fit(g = 2)
  expect_silent({
    V_naive <- vcov(fit, method = "naive")
    V_louis <- vcov(fit, method = "louis")
    V_oakes <- vcov(fit, method = "oakes")
    V_sand  <- vcov(fit, method = "sandwich")
  })
  for (V in list(V_naive, V_louis, V_oakes, V_sand)) {
    expect_true(is.matrix(V))
    expect_true(isSymmetric(V, tol = 1e-8))
    expect_true(all(is.finite(V)))
    eigs <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
    expect_true(min(eigs) > -tol_eig)
  }
})

test_that("vcov.csam works for g = 3 (multiple free pi parameters)", {
  fit <- make_and_fit(g = 3)
  expect_silent({
    V_naive <- vcov(fit, method = "naive")
    V_louis <- vcov(fit, method = "louis")
    V_oakes <- vcov(fit, method = "oakes")
    V_sand  <- vcov(fit, method = "sandwich")
  })
  for (V in list(V_naive, V_louis, V_oakes, V_sand)) {
    expect_true(is.matrix(V))
    expect_true(isSymmetric(V, tol = 1e-8))
    expect_true(all(is.finite(V)))
    eigs <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
    expect_true(min(eigs) > -tol_eig)
  }
})
