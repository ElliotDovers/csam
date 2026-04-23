// csam.cpp
#define TMBAD_FRAMEWORK
#include <TMB.hpp>
//#include "init.h" // for R CMD check: R_registerRoutines, R_useDynamicSymbols
#include "csam_utils.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(lik_type);

  // ---- Shared DATA ----
  DATA_MATRIX(Y);
  DATA_MATRIX(X);
  DATA_INTEGER(family);

  // ---- Shared PARAMETERS ----
  PARAMETER_VECTOR(beta0);
  PARAMETER_MATRIX(B);
  PARAMETER_VECTOR(logit_pi);
  PARAMETER_MATRIX(U);
  PARAMETER_MATRIX(Lambda);
  PARAMETER_VECTOR(log_phi);

  // ---- Dimensions ----
  int n = Y.rows();
  int s = Y.cols();
  int p = X.cols();
  int g = B.rows();
  int d = U.cols();

  // --- PARAMETER CONSTRAINTS

  // dispersion > 0
  vector<Type> phi = log_phi.exp();

  // 0 <= pi <= 1 (via softmax)
  // vector<Type> logpi(g);
  // {
  //   // subtract max for numerical stability
  //   Type m = logit_pi.maxCoeff();
  //
  //   // compute exp(theta - m)
  //   vector<Type> exp_theta = (logit_pi.array() - m).exp();
  //
  //   // compute log(sum(exp(theta - m)))
  //   Type log_denom = log(exp_theta.sum()) + m;
  //
  //   // final log(pi_k)
  //   for(int k = 0; k < g; k++) {
  //     logpi(k) = logit_pi(k) - log_denom;
  //   }
  // }

  // 0 <= pi <= 1 (via simplex w. stick-breaking)
  vector<Type> logpi(g);
  {
    int m = g - 1;

    vector<Type> logv(m);
    vector<Type> log1mv(m);

    // stable log(v_k) and log(1 - v_k)
    for(int k = 0; k < m; k++) {
      Type z = logit_pi(k);
      logv(k)   = -log1p(exp(-z));
      log1mv(k) = -log1p(exp(z));
    }

    // compute log(pi_k)
    Type acc = Type(0);
    for(int k = 0; k < m; k++) {
      logpi(k) = logv(k) + acc;
      acc += log1mv(k);
    }

    // final component
    logpi(g - 1) = acc;
  }
  vector<Type> pi = logpi.exp();
  REPORT(pi);
  // // Constraints on factor loadings
  // int idx = 0;
  // matrix<Type> Lambda_con(s, d);
  // Lambda_con.setZero();
  //
  // // Fill lower-triangular: for r = 0..d-1, allow entries j=r..s-1 (or choose per your convention)
  // for (int r = 0; r < d; ++r) {
  //   for (int j = r; j < s; ++j) {
  //     if (j == 0 && r == d) {
  //       // upper triangular:
  //       Lambda_con(j, r) = 0;
  //     } else if (j == r) {
  //       // diagonal: parameterise on log-scale to enforce positivity
  //       Lambda_con(j, r) = exp(Lambda(idx));
  //     } else {
  //       // off-diagonal: unconstrained
  //       Lambda_con(j, r) = Lambda(idx);
  //     }
  //     idx++;
  //   }
  // }
  matrix<Type> Lambda_con = Lambda;

  // ---- Penalised-only DATA ----
  DATA_SCALAR(psi1);
  DATA_SCALAR(psi2);

  // ---- Dispatch ----
  switch(lik_type) {

  case 0:
    #include "pll.h"
    break;

  case 1:
    #include "mll.h"
    break;

  default:
    error("lik_type must be 0 (penalised) or 1 (random effects).");
  }

  // This is never reached, but required syntactically
  return Type(0);
}
