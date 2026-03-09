#ifndef CSAM_UTILS_H
#define CSAM_UTILS_H

template<class Type>
Type log_sum_exp(const vector<Type> &v) {
  Type m = v.maxCoeff();
  return m + log((v.array() - m).exp().sum());
}

template<class Type>
Type linkinv(Type eta, int family) {
  switch(family) {
  case 0: return exp(eta);
  case 1: return Type(1)/(Type(1)+exp(-eta));
  case 2: return eta;
  case 3: return exp(eta);
  case 4: return exp(eta);
  default: error("Unknown family in linkinv");
  }
}

template<class Type>
Type loglik_y(Type y, Type mu, int family, Type phi=1.0) {
  switch(family) {
  case 0: return dpois(y, mu, true);
  case 1: return dbinom(y, Type(1), mu, true);
  case 2: return dnorm(y, mu, sqrt(phi), true);
  case 3: return dgamma(y, Type(1)/phi, Type(1)/(mu*phi), true);
  case 4: {
    Type var = mu + phi * mu * mu;
    return dnbinom2(y, mu, var, true);
  }
  default: error("Unknown family in loglik_y");
  }
}

#endif
