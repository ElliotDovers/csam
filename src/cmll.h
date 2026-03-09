// mll.h
{

  // ---- Complete data Parameters ----
  PARAMETER_MATRIX(Z);

  matrix<Type> XB(n, g);
  XB.setZero();
  for(int k=0;k<g;k++)
    for(int pp=0;pp<p;pp++)
      XB.col(k) += X.col(pp) * B(k,pp);

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
      loglik_k(k) = Z(j, k) * (log(pi(k)) + ll_k);
    }
    nll -= loglik_k.sum();
  }

  // U ~ N(0,1)
  for(int i=0;i<n;i++)
    for(int r=0;r<d;r++)
      nll -= dnorm(U(i,r), Type(0), Type(1), true);

  return nll;
}
