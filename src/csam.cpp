// csam.cpp
#include <TMB.hpp>
#include "init.h" // for R CMD check: R_registerRoutines, R_useDynamicSymbols
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
  PARAMETER_VECTOR(theta_pi);
  PARAMETER_MATRIX(U);
  PARAMETER_MATRIX(Lambda);
  PARAMETER_VECTOR(logphi);

  // ---- Dimensions ----
  int n = Y.rows();
  int s = Y.cols();
  int p = X.cols();
  int g = B.rows();
  int d = U.cols();

  // --- PARAMETER CONSTRAINTS

  // dispersion > 0
  vector<Type> phi = logphi.exp();

  // 0 <= pi <= 1 (via softmax)
  vector<Type> pi(g);
  {
    Type m = theta_pi.maxCoeff();
    vector<Type> exp_theta = (theta_pi.array() - m).exp();
    pi = exp_theta / exp_theta.sum();
  }

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

  case 2:
    #include "cpll.h"
    break;

  case 3:
    #include "cmll.h"
    break;

  default:
    error("lik_type must be 0 (penalised) or 1 (random effects).");
  }

  // This is never reached, but required syntactically
  return Type(0);
}
