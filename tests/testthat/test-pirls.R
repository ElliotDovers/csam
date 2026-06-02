# tests/testthat/test-pirls.R
library(testthat)

sims = 100
tol = 1e-6

test_that("P-IRLS scheme in C++ matches unpenalised (Binomial) GLM fit to at most tolerance", {
  params_glm.fit <- matrix(NA, nrow = sims, ncol = 3)
  params_pirls.cpp <- matrix(NA, nrow = sims, ncol = 3)
  for (seed in 1:sims) {
    set.seed(seed)
    # simulate data
    n=50;p=2
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    eta_true <- 0.5 + 1.2 * x1 - 0.8 * x2
    prob <- 1 / (1 + exp(-eta_true))
    y <- rbinom(n, size = 1, prob = prob)
    X = as.matrix(cbind(int = rep(1, length(x1)), x1, x2))


    ma <- glm.fit(X, y2, family = binomial())
    # mb <- glm2::glm.fit2(X, y2, family = binomial())
    mc <- penalised_glm_fit_C(X = X, y = y2, family = binomial(), lambda = 0, penalty_weights = c(0, 1, 1), tol = 1e-6)

    params_glm.fit[seed, 1:3] <- ma$coefficients
    params_pirls.cpp[seed, 1:3] <- mc$coefficients
  }

  # test
  expect_true(max(params_glm.fit - params_pirls.cpp) <= tol)
})

test_that("P-IRLS scheme in C++ matches unpenalised (Poisson) GLM fit to at most tolerance", {
  params_glm.fit <- matrix(NA, nrow = sims, ncol = 3)
  params_pirls.cpp <- matrix(NA, nrow = sims, ncol = 3)
  for (seed in 1:sims) {
    set.seed(seed)
    # simulate data
    n=50;p=2
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    eta_true <- 0.5 + 1.2 * x1 - 0.8 * x2
    y = rpois(n, exp(eta_true))
    X = as.matrix(cbind(int = rep(1, length(x1)), x1, x2))


    ma <- glm.fit(X, y, family = poisson())
    # mb <- glm2::glm.fit2(X, y2, family = poisson())
    mc <- penalised_glm_fit_C(X = X, y = y, family = poisson(), lambda = 0, penalty_weights = c(0, 1, 1), tol = 1e-6)

    params_glm.fit[seed, 1:3] <- ma$coefficients
    params_pirls.cpp[seed, 1:3] <- mc$coefficients
  }

  # test
  expect_true(max(params_glm.fit - params_pirls.cpp) <= tol)
})


