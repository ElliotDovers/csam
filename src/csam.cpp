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
  PARAMETER_VECTOR(phi);

  // ---- Penalised-only DATA ----
  DATA_SCALAR(psi1);
  DATA_SCALAR(psi2);

  // ---- Dimensions ----
  int n = Y.rows();
  int s = Y.cols();
  int p = X.cols();
  int g = B.rows();
  int d = U.cols();

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
