
library(TMB)
library(Rcpp)

tmb_binomial<- "
  #include <TMB.hpp>
  template<class Type>
  Type objective_function<Type>::operator () (){
    //data:
    DATA_MATRIX(Y);
    DATA_MATRIX(X);
    //parameters:
    PARAMETER_MATRIX(beta);
    // pre-calcs
    // rowsums of Y = n
    int p = Y.cols();
    matrix<Type> ones(p,1);
    ones << 1,1;
    matrix<Type> n = Y * ones; // rowsums of Y
    //procedure:
    int m = Y.rows();
    matrix<Type> eta = X * beta;
    matrix<Type> P = eta; 
    for(int i = 0; i < m; i++){ // C++ starts loops at 0!
        P(i,0) = exp(eta(i,0)) / (1.0 + exp(eta(i,0)));
    }
    Type nll = 0.0; // initialize negative log likelihood
    for(int i = 0; i < m; i++){ // C++ starts loops at 0!
      nll -= dbinom(Y(i,0), n(i,0), P(i,0), true);
      //nll -= dmultinom(Y(i), P(i), true);
    }
   return nll;
  }
"

write(tmb_binomial, file = "binomial_model.cpp")

compile("binomial_model.cpp")

dyn.unload(dynlib("binomial_model")) ## before loading to be safe
dyn.load(dynlib("binomial_model"))

Y <- t(rmultinom(10, size = 20, prob = c(0.2, 0.8)))
library(nnet)
mnom.fit <- multinom(Y~1)
X <- matrix(1, nrow = nrow(Y)) ##model.matrix(mnom.fit)

## strangely need to remove obj here if changing the data
rm(obj)

obj <- MakeADFun(
         data = list(Y = Y, X = X),
         parameters = list(beta = matrix(0)),
         DLL = "binomial_model")

(opt<-do.call("optim", obj))

plogis(opt$par)
predict(mnom.fit, type = "p")

