// pll.h
{
  // Softmax
  vector<Type> pi(g);
  {
    Type m = theta_pi.maxCoeff();
    vector<Type> exp_theta = (theta_pi.array() - m).exp();
    pi = exp_theta / exp_theta.sum();
  }

  // XB
  matrix<Type> XB(n, g);
  XB.setZero();
  for(int k=0;k<g;k++)
    for(int pp=0;pp<p;pp++)
      XB.col(k) += X.col(pp) * B(k,pp);

  // UL
  matrix<Type> UL(n, s);
  UL.setZero();
  for(int j=0;j<s;j++)
    for(int r=0;r<d;r++)
      UL.col(j) += U.col(r) * Lambda_con(j,r);

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

  // penalties
  nll += 0.5 * psi1 * (U.array()*U.array()).sum();
  nll += 0.5 * psi2 * (Lambda.array()*Lambda_con.array()).sum();

  return nll;
}

