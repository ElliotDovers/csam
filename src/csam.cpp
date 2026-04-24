// csam.cpp
#include <TMB.hpp>
#include <R_ext/Error.h>
#include "init.h"
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
  const int n = Y.rows();
  const int s = Y.cols();
  const int p = X.cols();
  const int g = B.rows();
  const int d = U.cols();

  // ---- Parameter constraints ----

  // dispersion > 0
  vector<Type> phi = exp(log_phi);

  // 0 <= pi_k <= 1 via stick-breaking simplex
  vector<Type> logpi(g);
  {
    const int m = g - 1;

    vector<Type> logv(m);
    vector<Type> log1mv(m);

    for (int k = 0; k < m; k++) {
      Type z = logit_pi(k);
      logv(k)   = -log1p(exp(-z));  // log(sigmoid(z))
      log1mv(k) = -log1p(exp(z));   // log(1 - sigmoid(z))
    }

    Type acc = Type(0);
    for (int k = 0; k < m; k++) {
      logpi(k) = logv(k) + acc;
      acc     += log1mv(k);
    }

    logpi(g - 1) = acc;
  }

  vector<Type> pi = exp(logpi);
  REPORT(pi);

  // ---- Factor loadings (currently unconstrained passthrough) ----
  matrix<Type> Lambda_con = Lambda;

  // ---- Penalised-only DATA ----
  DATA_SCALAR(psi1);
  DATA_SCALAR(psi2);

  // ---- Dispatch likelihood ----
  switch (lik_type) {

  case 0:
#include "pll.h"
    break;

  case 1:
#include "mll.h"
    break;

  default:
    error("lik_type must be 0 (penalised) or 1 (random effects).");
  }

  // Required syntactically
  return Type(0);
}
