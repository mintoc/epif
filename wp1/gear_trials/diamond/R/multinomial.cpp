#include <TMB.hpp>                                

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Y);                                 // Response counts
  DATA_MATRIX(X);                                 // Model matrix
  PARAMETER_VECTOR(theta);                        // pars
  // get the number of groups
  int m = Y.cols();
  std::cout << "Number of groups = " << m << "\n";
  // get the number of parameters per group (except first)
  int p = X.cols();
  std::cout << "Number of parameters per group = " << p << "\n";
  // create the parameter matrix
  //matrix<double> beta(p,m); 
  matrix<Type> beta(p,m); 
  // leave first column at zero
  beta.block(0,1,2,3) << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::cout << "Dimensions beta \n" << beta << "\n";

  Type ans;
  return ans;
}
