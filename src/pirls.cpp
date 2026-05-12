#include <RcppEigen.h>
using Eigen::MatrixXd;
using Eigen::VectorXd;

inline double binomial_deviance(const VectorXd& y, const VectorXd& mu) {
  return -2.0 * (
      y.array() * mu.array().log() +
        (1.0 - y.array()) * (1.0 - mu.array()).log()
  ).sum();
}

inline double poisson_deviance(const VectorXd& y, const VectorXd& mu) {
  double dev = 0.0;
  for (int i = 0; i < y.size(); ++i) {
    if (y[i] > 0.0)
      dev += y[i] * std::log(y[i] / mu[i]) - (y[i] - mu[i]);
    else
      dev += mu[i];
  }
  return 2.0 * dev;
}

extern "C" SEXP penalised_glm_irls_cpp(
    SEXP X_, SEXP y_, SEXP weights_, SEXP offset_,
    SEXP penalty_weights_, SEXP lambda_,
    SEXP beta_init_, SEXP family_,
    SEXP maxit_, SEXP tol_
) {
  BEGIN_RCPP;

  const MatrixXd X = Rcpp::as<MatrixXd>(X_);
  const VectorXd y = Rcpp::as<VectorXd>(y_);
  const VectorXd w = Rcpp::as<VectorXd>(weights_);
  const VectorXd offset = Rcpp::as<VectorXd>(offset_);
  const VectorXd pen_w = Rcpp::as<VectorXd>(penalty_weights_);
  const double lambda = Rcpp::as<double>(lambda_);
  VectorXd beta = Rcpp::as<VectorXd>(beta_init_);

  const int family = Rcpp::as<int>(family_);
  const int maxit = Rcpp::as<int>(maxit_);
  const double tol = Rcpp::as<double>(tol_);

  const int n = X.rows();
  //const int p = X.cols();

  // ---- Constants ----
  constexpr double ETA_MAX = 30.0; // Poisson log-link safety bound
  constexpr double ETA_MIN = -30.0; // Poisson log-link safety bound
  //constexpr double MU_EPS = 1e-12; // Binomial safety bound to clip mu
  constexpr double MU_ETA_FLOOR = 1e-12;

  VectorXd eta(n), mu(n), mu_eta(n), var_mu(n);
  VectorXd W(n), z(n);

  bool converged = false;
  int iter = 0;

  for (iter = 0; iter < maxit; ++iter) {

    // ---- Current linear predictor ----
    eta = X * beta + offset;

    // ---- Mean & deviance at current beta ----
    double dev_old;

    if (family == 0) { // Poisson
      eta = eta.array().min(ETA_MAX).max(ETA_MIN);
      mu = eta.array().exp();
      dev_old = poisson_deviance(y, mu);
      mu_eta = mu;
      var_mu = mu;
    } else if (family == 1) { // Binomial
      mu = 1.0 / (1.0 + (-eta.array()).exp());
      //mu = mu.array().max(MU_EPS).min(1.0 - MU_EPS); // numerical clip causes problems in the EM algorithm
      dev_old = binomial_deviance(y, mu);
      mu_eta = mu.array() * (1.0 - mu.array());
      var_mu = mu_eta;
    } else {
      Rcpp::stop("Family not implemented");
    }

    mu_eta = mu_eta.array().max(MU_ETA_FLOOR);

    // ---- IRLS weights and working response ----
    W = w.array() * mu_eta.array().square() / var_mu.array();
    z = (eta - offset).array() + (y - mu).array() / mu_eta.array();

    // ---- Penalised normal equations ----
    MatrixXd XtWX = X.transpose() * (W.asDiagonal() * X);
    XtWX.diagonal().array() += lambda * pen_w.array();
    VectorXd XtWz = X.transpose() * (W.array() * z.array()).matrix();

    Eigen::LLT<MatrixXd> llt(XtWX);
    if (llt.info() != Eigen::Success) break;

    VectorXd beta_full = llt.solve(XtWz);

    // ---- Deviance-based line search ----
    double dev_new = dev_old;  // initialise defensively
    double step = 1.0;
    VectorXd beta_new = beta;
    bool accepted = false;

    while (step > 1e-6) {

      VectorXd beta_try = beta + step * (beta_full - beta);
      eta = X * beta_try + offset;

      if (family == 0) {
        // --- Poisson ---
        eta = eta.array().min(ETA_MAX).max(ETA_MIN);
        mu = eta.array().exp();

        if (!mu.array().isFinite().all()) {
          step *= 0.5;
          continue;
        }

        dev_new = poisson_deviance(y, mu);

      } else if (family == 1) {
        // --- Binomial ---
        mu = 1.0 / (1.0 + (-eta.array()).exp());
        //mu = mu.array().max(MU_EPS).min(1.0 - MU_EPS); // numerical clip causes problems in the EM algorithm

        if (!mu.array().isFinite().all()) {
          step *= 0.5;
          continue;
        }

        dev_new = binomial_deviance(y, mu);

      } else {
        Rcpp::stop("Family not implemented in pirls.cpp");
      }

      // ---- Accept only if deviance decreases ----
      if (dev_new <= dev_old) {
        beta_new = beta_try;
        accepted = true;
        break;
      }

      step *= 0.5;
    }

    if (!accepted) {
      // No acceptable step found: abort IRLS
      break;
    }

    // ---- Convergence check ----
    if (dev_old - dev_new < tol) {
      beta = beta_new;
      converged = true;
      break;
    }

    beta = beta_new;
  }

  return Rcpp::List::create(
    Rcpp::Named("coefficients") = beta,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations")  = iter + 1
  );

  END_RCPP;
}
