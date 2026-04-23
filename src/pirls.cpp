#include <RcppEigen.h>

// We rely on RcppEigen only (NO TMB, NO Rcpp.h)
using Eigen::MatrixXd;
using Eigen::VectorXd;

extern "C" SEXP penalised_glm_irls_cpp(
    SEXP X_,               // numeric matrix (n x p)
    SEXP y_,               // numeric vector (n)
    SEXP weights_,         // numeric vector (n)
    SEXP offset_,          // numeric vector (n)
    SEXP penalty_weights_, // numeric vector (p)
    SEXP lambda_,          // scalar
    SEXP beta_init_,       // numeric vector (p)
    SEXP family_,          // integer family code
    SEXP maxit_,           // integer
    SEXP tol_              // numeric
) {
  BEGIN_RCPP;

  // ---- Map R objects to Eigen ----
  const MatrixXd X      = Rcpp::as<MatrixXd>(X_);
  const VectorXd y      = Rcpp::as<VectorXd>(y_);
  const VectorXd w      = Rcpp::as<VectorXd>(weights_);
  const VectorXd offset = Rcpp::as<VectorXd>(offset_);
  const VectorXd pen_w  = Rcpp::as<VectorXd>(penalty_weights_);

  const double lambda   = Rcpp::as<double>(lambda_);
  VectorXd beta         = Rcpp::as<VectorXd>(beta_init_);
  const int family      = Rcpp::as<int>(family_);
  const int maxit       = Rcpp::as<int>(maxit_);
  const double tol      = Rcpp::as<double>(tol_);

  const int n = X.rows();
  const int p = X.cols();

  // ---- Constants ----
  constexpr double ETA_MAX = 30.0;     // Poisson log-link safety bound
  constexpr double ETA_MIN = -30.0;
  constexpr double MU_ETA_FLOOR = 1e-12;
  constexpr double RIDGE_FLOOR = 0.0;

  // ---- Working vectors ----
  VectorXd eta(n), mu(n), mu_eta(n), var_mu(n);
  VectorXd W(n), z(n);

  bool converged = false;
  int iter = 0;

  // ---- IRLS loop ----
  for (iter = 0; iter < maxit; ++iter) {

    // Save previous iterate
    VectorXd beta_old = beta;

    // ---- Linear predictor at current beta ----
    eta = X * beta_old + offset;

    // ---- Family-specific mean / variance ----
    if (family == 0) {
      // Poisson (log link)
      eta = eta.array().min(ETA_MAX).max(ETA_MIN);
      mu     = eta.array().exp();
      mu_eta = mu;
      var_mu = mu;

    } else if (family == 1) {
      // Binomial (logit link, size = 1)
      mu     = 1.0 / (1.0 + (-eta.array()).exp());
      mu_eta = mu.array() * (1.0 - mu.array());
      var_mu = mu_eta;

    } else {
      Rcpp::stop("Family not implemented in penalised_glm_irls_cpp");
    }

    // ---- Numerical guards ----
    mu_eta = mu_eta.array().max(MU_ETA_FLOOR);

    // ---- IRLS weights ----
    W = w.array() * mu_eta.array().square() / var_mu.array();

    // ---- Working response ----
    z = ( (eta - offset).array()
            + (y - mu).array() / mu_eta.array()
    ).matrix();


    // ---- Weighted crossproducts ----
    MatrixXd WX = X;
    WX.array().colwise() *= W.array();

    MatrixXd XtWX = X.transpose() * WX;
    VectorXd XtWz = X.transpose() * (W.array() * z.array()).matrix();

    // ---- Ridge penalty (+ small numerical floor) ----
    XtWX.diagonal().array() += lambda * pen_w.array();
    XtWX.diagonal().array() += RIDGE_FLOOR;

    // ---- Solve normal equations ----
    Eigen::LLT<MatrixXd> llt(XtWX);
    if (llt.info() != Eigen::Success) {
      break;  // numerical failure
    }

    VectorXd beta_full = llt.solve(XtWz);

    // ---- Try full Newton/IRLS step FIRST ----
    bool accept_full = beta_full.array().isFinite().all();

    if (accept_full) {
      eta = X * beta_full + offset;

      if (family == 0) {
        // Poisson: reject if eta explodes
        if (eta.array().abs().maxCoeff() > ETA_MAX) {
          accept_full = false;
        } else {
          mu = eta.array().exp();
          if (!mu.array().isFinite().all()) {
            accept_full = false;
          }
        }
      } else if (family == 1) {
        mu = 1.0 / (1.0 + (-eta.array()).exp());
        if (!mu.array().isFinite().all()) {
          accept_full = false;
        }
      }
    }

    // ---- Accept full step if safe ----
    if (accept_full) {
      beta = beta_full;

    } else {
      // ---- Step-halving fallback (only when needed) ----
      double step = 1.0;
      VectorXd beta_try(p);

      while (step > 1e-6) {
        beta_try = beta_old + step * (beta_full - beta_old);

        eta = X * beta_try + offset;

        if (family == 0) {
          eta = eta.array().min(ETA_MAX).max(ETA_MIN);
          mu  = eta.array().exp();
        } else if (family == 1) {
          mu = 1.0 / (1.0 + (-eta.array()).exp());
        }

        if (beta_try.array().isFinite().all() &&
            mu.array().isFinite().all()) {
          beta = beta_try;
          break;
        }

        step *= 0.5;
      }
    }

    // ---- Convergence check (correct reference) ----
    if ((beta - beta_old).cwiseAbs().maxCoeff() < tol) {
      converged = true;
      break;
    }
  }

  // ---- Return results ----
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = beta,
    Rcpp::Named("converged")    = converged,
    Rcpp::Named("iterations")  = iter + 1
  );

  END_RCPP;
}

