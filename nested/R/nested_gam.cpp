#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // FIRST LEVEL
  DATA_MATRIX(Y1);
  DATA_MATRIX(X1);
  DATA_MATRIX(Z1); // basis function
  DATA_MATRIX(O1); // offset
  DATA_MATRIX(G1); // level ID
  DATA_MATRIX(X1pred);
  DATA_MATRIX(Z1pred);
  // SECOND LEVEL
  DATA_MATRIX(Y2);
  DATA_MATRIX(X2);
  DATA_MATRIX(Z2); // basis function
  DATA_MATRIX(O2); // offset
  DATA_MATRIX(G2); // level ID  
  DATA_MATRIX(X2pred);
  DATA_MATRIX(Z2pred);  
  // pars
  // FIRST LEVEL
  PARAMETER_MATRIX(beta1);
  PARAMETER(lnsigmau1);
  PARAMETER_MATRIX(u1);
  PARAMETER(lnsigmab1);
  PARAMETER_MATRIX(b1);
  // SECOND LEVEL
  PARAMETER_MATRIX(beta2);
  PARAMETER(lnsigmau2);
  PARAMETER_MATRIX(u2);
  PARAMETER(lnsigmab2);
  PARAMETER_MATRIX(b2);  
  //---------------
  // LEVEL 1 CALCS
  //---------------
  // prelim calcs
  vector<Type> n = Y1.rowwise().sum();  
  Type sigmau1 = exp(lnsigmau1);
  Type sigmab1 = exp(lnsigmab1);
  Type sigmab2 = exp(lnsigmab2);
  Type nll = 0;
  // spline component
  // second order penalty
  int m1 = Z1.cols();
  for(int i = 2; i < m1; i++){
    nll -= dnorm((u1(i,0) - u1(i-1,0)), (u1(i-1,0) - u1((i-2),0)), sigmau1, true);
  }
  // center the splines
  nll -= dnorm(sum(u1), Type(0.0), Type(1.0), true);
  // haul-level 1 random effects
  nll -= sum(dnorm(vector<Type>(b1), Type(0.0), sigmab1, true));
  // observations
  matrix<Type> f1 = Z1 * u1;
  //f1 /= sigmau1;
  vector<Type> eta1 = X1 * beta1 + f1 + G1 * b1 + O1;
  vector<Type> p1 = invlogit(eta1);
  // NB n which is the total count used here
  nll -= sum(dbinom(vector<Type>(Y1.col(0)), n, p1, true));
  //---------------
  // LEVEL 2 CALCS
  //---------------
  // prelim calcs
  Type sigmau2 = exp(lnsigmau2);
  // spline component
  // second order penalty
  int m2 = Z2.cols();
  for(int i = 2; i < m2; i++){
    //nll -= dnorm((u2(i,0) - 2.0 * u2(i-1,0) + u2((i-2),0)), Type(0.0), sigmau2, true);
    nll -= dnorm((u2(i,0) - u2(i-1,0)), (u2(i-1,0) - u2((i-2),0)), sigmau2, true);
  }
  // center the splines
  nll -= dnorm(sum(u2), Type(0.0), Type(1.0), true);
  // haul-level 2 random effects
  nll -= sum(dnorm(vector<Type>(b2), Type(0.0), sigmab2, true));  
  // observations
  matrix<Type> f2 = Z2 * u2;
  //f2 /= sigmau2;
  vector<Type> eta2 = X2 * beta2 + f2 + G2 * b2 + O2;
  vector<Type> p2 = invlogit(eta2);
  vector<Type> n2 = Y2.rowwise().sum();  
  // NB n which is the total count used here
  vector<Type> pL = p1.array() * p2.array();
  vector<Type> etaL = logit(pL);
  vector<Type> pU = p1.array() * (Type(1.0) - p2.array());
  vector<Type> etaU = logit(pU);
  nll -= sum(dbinom(vector<Type>(Y2.col(0)), n, pL, true));
  // predictions
  // P(T)
  matrix<Type> f1pred = Z1pred * u1;
  //f1pred /= sigmau1;
  vector<Type> eta1pred = X1pred * beta1 + f1pred;  
  // P(L|T)
  matrix<Type> f2pred = Z2pred * u2;
  //f2pred /= sigmau2;
  vector<Type> eta2pred = X2pred * beta2 + f2pred;
  ADREPORT(eta1pred);
  ADREPORT(eta2pred);
  return nll;
}
