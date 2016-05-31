
## generate some linear random effects
## using rmvnorm 
## include

library(TMB)
library(Rcpp)

##----------
## DATA SIM 
##----------

## fixed effects
n <- 100 ## obs per group
m <- 10 ## groups

beta0 <- matrix(rnorm(3*2, sd = 0.1), ncol = 3, nrow = 2)
beta <- cbind(0, beta0)
##X <- cbind(1, rnorm(n))
X0 <- cbind(1, seq(-2,2, length = n))
X <- X0[rep(1:n, times = m), ]

## random effects
library(mvtnorm)
Sigma <- matrix(0, 3, 3)

Sigma[1,1] <- 0.07##0.2
Sigma[2,2] <- 0.07##0.4
Sigma[3,3] <- 0.07##0.1

u <- rmvnorm(m, mean = c(0, 0, 0), sigma = Sigma)
u <- cbind(1:m, u)
u.mat <- u[rep(1:m, each = n), ]

## two conditional variables
Xcond <- array(rnorm(nrow(X) * 4 * 2), dim = c(nrow(X), 4, 2))

beta.Xcond <- rnorm(2, sd = 0.1)

eta <- X %*% beta + cbind(0, u.mat[, -1]) + beta.Xcond[1] * Xcond[,,1] + beta.Xcond[2] * Xcond[,,2]
P <- exp(eta) / rowSums(exp(eta)) 

## Counts
Y <- t(sapply(1:nrow(P), function(z){rmultinom(1, size = 50, prob = P[z,])}))

P.df <- data.frame(gp = c(rep(u.mat[,1] - 1, times = 4)),
                   category = rep(1:4, each = nrow(P)),
                   x = rep(X[, 2], 4),
                   proportion = c(P),
                   obs.proportion = c(prop.table(Y, margin = 1)))

## library(ggplot2)
## theme_set(theme_bw())
## ggplot(data = P.df, aes(x = x, y = proportion, colour = factor(category))) + geom_line() + geom_point(aes(y = obs.proportion)) + facet_grid(gp ~ category) + theme(legend.position = "none")

##-------------
## TMB FITTING 
##-------------

tmb_multinomial_re<- "
  #include <TMB.hpp>
  using namespace density;
  //using namespace std;
  template <class Type>
  Type objective_function<Type>::operator () (){
    //------
    // DATA
    //------
    DATA_MATRIX(Y);
    DATA_MATRIX(X);
    DATA_ARRAY(Xcond);
    DATA_IVECTOR(gp);
    DATA_INTEGER(ngp);
    //------------
    // PARAMETERS
    //------------
    PARAMETER_MATRIX(beta0);
    PARAMETER_VECTOR(betacond);
    PARAMETER_VECTOR(a); // Cholesky factors
    PARAMETER_MATRIX(u0);
    //---------------------------
    // PRELIMINRARY CALCULATIONS
    //---------------------------
    // NUMBER OF OBSERVATIONS
    int n = Y.rows();
    // NUMBER OF CHOICES
    int m = Y.cols();
    // number of covariates
    int p = X.cols();
    // depth of the array
    int q = Xcond.cols(); // number of outermost dimensions - slabs
    // BETA WITH A COLUMN OF ZEROS FIRST
    matrix<Type> beta(p,m);
    for(int i = 0; i < p; i++){
      for(int j = 1; j < m; j++){ // leaves the first column as zeros
        beta(i,j) = beta0(i,j-1);
      }
    }
    // RANDOM EFFECTS WITH COLUMN OF ZEROS FIRST
    matrix<Type> u1(ngp,m);
    for(int i = 0; i < ngp; i++){
      for(int j = 1; j < m; j++){ // leaves the first column as zeros
        u1(i,j) = u0(i,j-1);
      }
    }
    // full matrix of random effects
    matrix<Type> u(n,m);
    for(int i = 0; i < n; i++){
        u.row(i) = u1.row(gp[i]);
    }
    // sigma
    matrix<Type> L(3,3);
    int ii = 0;
    for(int i = 0; i < 3; i++){
      for(int j = 0; j <= i; j++){
        L(i,j) = a(ii);
        ii += 1;
      }
    }
    matrix<Type> Sigma = L * L.transpose();
    //-----------
    // PROCEDURE
    //-----------
    matrix<Type> eta = X * beta + u;
    // conditional variables
    for(int i = 0; i < q; i++){
      matrix<Type> etacondi = asMatrix(vector<Type>(Xcond.col(i) * betacond(i)), m, n); // transposed below 
      //std::cout << etacondi.transpose() << std::endl << std::endl;
      eta += etacondi.transpose();
    }
    matrix<Type> expeta = eta;
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
          expeta(i,j) = exp(eta(i,j)); // can we do quicker?
      }
    }
    matrix<Type> ones(m,1);
    ones.fill(1.0);
    matrix<Type> rowsumsexpeta =  expeta * ones;
    matrix<Type> P = eta;
    for(int i = 0; i < n; i++){ 
      for(int j = 0; j < m; j++){ 
        P(i,j) = expeta(i,j) / rowsumsexpeta(i,0);
      }
    }
    Type nll = 0.0; // initialize negative log likelihood
    //vector<Type> yrow(m), prow(m);
    // random effects component
    for(int j = 0; j < ngp; j++){
      nll += MVNORM(Sigma, false)(vector<Type>(u0.row(j))); // note + as return negative log density
    }
    // observations component
    for(int i = 0; i < n; i++){ 
      // cast the rows to vectors and use dmultinom
      //;
      //;
      //yrow = Y.row(i);
      //prow = P.row(i);
      //nll -= dmultinom(yrow, prow, true);
      nll -= dmultinom(vector<Type>(Y.row(i)), vector<Type>(P.row(i)), true);
    }
   return nll;
  }
"

write(tmb_multinomial_re, file = "multinomial_re_model.cpp")

compile("multinomial_re_model.cpp")

dyn.unload(dynlib("multinomial_re_model")) ## before loading to be safe
dyn.load(dynlib("multinomial_re_model"))

## strangely need to remove obj here if changing the data
rm(obj)

gp <- u.mat[,1] - 1
ngp <- length(unique(gp))

obj <- MakeADFun(
    data = list(Y = Y, X = X, Xcond = Xcond, gp = gp, ngp = ngp),
    parameters = list(
        beta0 = matrix(0, ncol = ncol(Y) - 1, nrow = ncol(X)),
        betacond = rep(0, dim(Xcond)[3]),
        a = rep(0.1, 6),
        u0 = matrix(0, ncol = ncol(Y) - 1, nrow = ngp)
        ),
    random = c("u0"),
    DLL = "multinomial_re_model")

(opt<-do.call("optim", obj))

rep <- sdreport(obj)

fixed.names <- names(rep$par.fixed)

(beta.hat <- matrix(rep$par.fixed[grep("beta", fixed.names)], nrow = ncol(X)))

a.hat <- rep$par.fixed[fixed.names == "a"]
L.hat <- matrix(0, 3, 3)
ii <- 1
for(i in 1:3){
    for(j in 1:i){
        L.hat[i,j] <- a.hat[ii]
        ii <- ii + 1
    }
}

Sigma.hat <- L.hat %*% t(L.hat)

u.hat <- matrix(rep$par.random, ncol = 3)

plot(as.data.frame(u.hat))

plot(u.hat, u[, -1]); abline(c(0,1))

