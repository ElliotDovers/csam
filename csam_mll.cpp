// csam.cpp
#include <TMB.hpp>
//#include "init.h" // for R CMD check: R_registerRoutines, R_useDynamicSymbols
#include "csam_utils.h"

template<class Type>
  Type objective_function<Type>::operator() ()
{
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

  // 0 <= pi <= 1 & sum(pi) = 1 (via softmax)
  // vector<Type> pi(g);
  // {
  //   Type m = theta_pi.maxCoeff();
  //   vector<Type> exp_theta = (theta_pi.array() - m).exp();
  //   pi = exp_theta / exp_theta.sum();
  // }

  // 0 <= pi <= 1 & sum(pi) = 1 (via inverse logit and g-1 free parameters, so-call stick break or simplex)
  vector<Type> pi(g);
  {
    // v_k = invlogit(theta_pi_k) for k = 0..g-2
    vector<Type> v(g-1);
    for(int k=0; k < g-1; k++) {
      v(k) = Type(1) / (Type(1) + exp(-theta_pi(k)));
    }

    // stick breaking: pi_k = v_k * prod_{l<k} (1 - v_l)
    Type remaining = Type(1);
    for(int k=0; k < g-1; k++) {
      pi(k) = v(k) * remaining;
      remaining *= (Type(1) - v(k));
    }
    // last component gets the remaining stick
    pi(g-1) = remaining;
  }


  matrix<Type> XB(n, g);
  XB.setZero();
  for(int k=0;k<g;k++)
    for(int pp=0;pp<p;pp++)
      XB.col(k) += X.col(pp) * B(k,pp);

  matrix<Type> UL(n, s);
  UL.setZero();
  for(int j=0;j<s;j++)
    for(int r=0;r<d;r++)
      UL.col(j) += U.col(r) * Lambda(j,r);

  Type nll = 0.0;

  for(int j=0;j<s;j++) {
    vector<Type> loglik_k(g);
    loglik_k.setZero();

    for(int k=0;k<g;k++) {
      Type ll_k = 0.0;
      for(int i=0;i<n;i++) {
        Type eta = beta0(j) + XB(i,k) + UL(i,j);
        Type mu  = linkinv<Type>(eta, family);
        ll_k += loglik_y<Type>(Y(i,j), mu, family, phi(j));
      }
      loglik_k(k) = log(pi(k)) + ll_k;
    }
    nll -= log_sum_exp<Type>(loglik_k);
  }

  // U ~ N(0,1)
  for(int i=0;i<n;i++)
    for(int r=0;r<d;r++)
      nll -= dnorm(U(i,r), Type(0), Type(1), true);
  REPORT(pi);
  return nll;
}
