
library(TMB)
library(Rcpp)

tmb_multinomial<- "
  #include <TMB.hpp>
  template <class Type>
  Type objective_function<Type>::operator () (){
    //------
    // DATA
    //------
    DATA_MATRIX(Y);
    DATA_MATRIX(X);
    //------------
    // PARAMETERS
    //------------
    PARAMETER_MATRIX(beta0);
    //---------------------------
    // PRELIMINRARY CALCULATIONS
    //---------------------------
    // NUMBER OF OBSERVATIONS
    int n = Y.rows();
    // NUMBER OF CHOICES
    int m = Y.cols();
    // number of covariates
    int p = X.cols();
    // BETA WITH A COLUMN OF ZEROS FIRST
    matrix<Type> beta(p,m);
    for(int i = 0; i < p; i++){
      for(int j = 1; j < m; j++){ // leaves the first column as zeros
        beta(i,j) = beta0(i,j-1);
      }
    }
    //-----------
    // PROCEDURE
    //-----------
    matrix<Type> eta = X * beta;
    matrix<Type> expeta = eta;
    for(int i = 0; i < n; i++){ // C++ starts loops at 0!
      for(int j = 0; j < m; j++){ // C++ starts loops at 0!
          expeta(i,j) = exp(eta(i,j)); // can we do quicker?
      }
    }
    matrix<Type> ones(m,1);
    ones.fill(1.0);
    matrix<Type> rowsumsexpeta =  expeta * ones;
    matrix<Type> P = eta;
    for(int i = 0; i < n; i++){ // C++ starts loops at 0!
      for(int j = 0; j < m; j++){ // C++ starts loops at 0!
        P(i,j) = expeta(i,j) / rowsumsexpeta(i,0);
      }
    }
    Type nll = 0.0; // initialize negative log likelihood
    vector<Type> yrow(m), prow(m);
    for(int i = 0; i < n; i++){ // C++ starts loops at 0!
      // cast the rows to vectors and use dmultinom
      //vector<Type> yrow = Y.row(i);
      //vector<Type> prow = P.row(i);
      yrow = Y.row(i);
      prow = P.row(i);
      nll -= dmultinom(yrow, prow, true);
    }
   return nll;
  }
"

write(tmb_multinomial, file = "multinomial_model.cpp")

compile("multinomial_model.cpp")

dyn.unload(dynlib("multinomial_model")) ## before loading to be safe
dyn.load(dynlib("multinomial_model"))

n <- 10
Y <- t(rmultinom(n, size = 20, prob = c(0.3, 0.2, 0.4, 0.1)))
library(nnet)

X <- cbind(1, rnorm(n), rnorm(n))
mnom.fit <- multinom(Y ~ X[, -1])

## strangely need to remove obj here if changing the data
rm(obj)

obj <- MakeADFun(
         data = list(Y = Y, X = X),
         parameters = list(beta0 = matrix(0, ncol = ncol(Y) - 1, nrow = ncol(X))),
         DLL = "multinomial_model")

(opt<-do.call("optim", obj))

t(matrix(opt$par, ncol = ncol(Y) - 1, nrow = ncol(X)))
coef(mnom.fit)

Phat <- predict(mnom.fit, type = "prob")

sum(sapply(1:nrow(Y), function(z){dmultinom(Y[z,], size = 20, prob = Phat[z,], log = TRUE)}))

