library(nnet)
library(TMB)

setwd("C:/Users/bburke/Desktop/GMIT/Gear Trials/R/multinomial_re_ADMB_TMB/newdir")

##Download data (http://www.nuigalway.ie/faculties_departments/mathematics/staff/John_Hinde/SMIR/SMIR.html)
##and load it in 
load("C:/Users/bburke/Desktop/GMIT/R/Statistics in R - J Hinde/SMIR/data/R_data/miners.rda")

mnomfit <- multinom(count ~ period + resp + lyr + resp:lyr, data=miners.long)

Y <- cbind(miners$n, miners$m, miners$s)
respm.lyr <- cbind(rbind(rep(0,8),miners.long$lyr[1:8],rep(0,8)))
resps.lyr <- cbind(rbind(rep(0,8),rep(0,8),miners.long$lyr[1:8]))
X <- cbind(miners.long$period,miners.long$resp,miners.long$lyr,resps.lyr,respm.lyr)
##X <- model.matrix(miners.glm8)
X <- model.matrix(mnomfit)
rownames(X)<- NULL

##-----
## TMB 
##-----
## Model excludes intercept
## model spec
tmb_multinomial_cs<-"
#include <TMB.hpp>
using namespace density;
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
  //PARAMETER(alpha);
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
//  for(int i = 0; i < n; i++){
//    for(int j = 0; j < m; j++){
//        eta(i,j) = eta(i,j) + alpha; // can we do quicker?
//      }
//    }
// element-wise exponentiation
  matrix<Type> expeta(n,m);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
        expeta(i,j) = exp(eta(i,j)); // can we do quicker?
      }
    }
// row sums of exp(eta)
  matrix<Type> ones(m,1);
  ones.fill(1.0);
  matrix<Type> rowsumsexpeta =  expeta * ones;

  matrix<Type> P(n,m);
    for(int i = 0; i < n; i++){ 
      for(int j = 0; j < m; j++){ 
        P(i,j) = expeta(i,j) / rowsumsexpeta(i,0);
      }
    }

Type nll = 0.0; // initialize negative log likelihood

// observations component
  for(int i = 0; i < n; i++){ 
  nll -= dmultinom(vector<Type>(Y.row(i)), vector<Type>(P.row(i)), true);
  }
return nll;
}
"

write(tmb_multinomial_cs, file = "multinomial_cs.cpp")

compile("multinomial_cs.cpp")

dyn.unload(dynlib("multinomial_cs")) ## before loading to be safe
dyn.load(dynlib("multinomial_cs"))

rm(obj)

obj <- MakeADFun(
  data = list(
    Y = Y,
    X = X),
  parameters = list(
    ##alpha = 0,
    beta0 = matrix(0, ncol = ncol(Y) - 1, nrow = ncol(X))
    ),
  map = list(
    ## if need to fix any parameters
    ##betacond  = factor(rep(NA, dim(Xcond.array)[3])),
    ##gammaw = factor(rep(NA, 4)),
    ##gammawcl = factor(rep(NA, 4))
    ##a = factor(rep(NA, 6)),
    ##u0 = factor(matrix(NA, ncol = ncol(Y) - 1, nrow = ngp))
  ),
  DLL = "multinomial_cs"
)

(opt<-do.call("optim", obj))

(rep <- sdreport(obj))

##Relevel the resp variable
miners.long$respm <- relevel(miners.long$resp,ref="m")

##Reshape/Create variable that compresses three categories into two
miners.long <- transform(miners.long,r23=ifelse(resp!="n",1,0))
miners.long <- transform(miners.long,r23lyr=r23*lyr)

##plot
plot(opt$par, xlim=c(0,11), ylim=c(-11,10))
points(coef(miners.glm8), col="blue")
