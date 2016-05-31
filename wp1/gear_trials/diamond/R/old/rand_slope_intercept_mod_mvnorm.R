##setwd("C:/Users/bburke/Desktop/GMIT/Gear Trials/R/TMB")
library(TMB)
library(Rcpp)
library(lme4)

## moved data part here
m <- 20
n <- 10
alpha<- 2
beta <- 1
sigma.b <- 3 ## random effects slope
sigma.n <- 1 ## residual SD
sigma.a <- 3 ## random effects intercept

obs <- m * n

a <- rnorm(m, mean = alpha, sd=sigma.a)
b <- rnorm(m, mean = beta, sd = sigma.b)

y.mat <- matrix(NA, nrow = n, ncol = m)
x<- rnorm(n)

for(i in 1:m){
  y.mat[, i] <- rnorm(n, mean =a[i]+b[i]*x, sd = sigma.n)
}

X <- rep(x, times=m)

df <- data.frame(gps = rep(0:(m-1),each = n), y = c(y.mat), x=X)
df$fgps<-as.factor(df$gps)

##library(ggplot2)
##ggplot(df, aes(x=x, y=y)) + geom_point() + facet_wrap(~gps)

Y <- matrix(df$y, ncol = 1)

tmb_rand<-"
#include<TMB.hpp>
using namespace density;
template<class Type>
Type objective_function<Type>::operator()(){
//DATA:
  DATA_VECTOR(Y);
  DATA_VECTOR(X);
  DATA_IVECTOR(gp);
  DATA_SCALAR(m);//groups
  DATA_SCALAR(n);//individuals
//PARAMETER:
  PARAMETER(alpha);
  PARAMETER(beta);
  PARAMETER_VECTOR(cholfactor); // optimize on the real line and then transform
  PARAMETER(logsigma_n);
//PROCEDURE:
  PARAMETER_MATRIX(theta); // random effects matrix
  matrix<Type> L(2,2);
  L(0,0) = cholfactor[0];
  L(1,0) = cholfactor[1];
  L(1,1) = cholfactor[2];
  matrix<Type> Sigma = L * L.transpose();
  Type sigma_n = exp(logsigma_n);
//
 //std::cout << MVNORM(Sigma, false)(test1) << std::endl << std::endl;
Type nll = 0.0;
for(int j = 0; j < m; j++){
  nll += MVNORM(Sigma, false)(vector<Type>(theta.row(j)));
}
//
for(int i = 0; i < n; i++){
  int indexi = gp[i];
  nll -= dnorm(Y[i], (alpha + theta(indexi,0)) + (beta + theta(indexi,1)) * X[i], sigma_n, true);
}
return nll;
}
"

write(tmb_rand, file = "rand_model.cpp")

compile("rand_model.cpp")

dyn.unload(dynlib("rand_model")) ## before loading to be safe
dyn.load(dynlib("rand_model"))

rm(obj)

## note random statement

obj <- MakeADFun(
  data = list(Y = df$y, X = df$x , n = nrow(Y), m = m, gp = df$gps),
  parameters = list(alpha = 0, beta = 0, cholfactor = c(1, 0, 1), logsigma_n = 0, theta = matrix(0, nrow = m, ncol = 2)),
  random = c("theta"),
  DLL = "rand_model")

(opt<-do.call("optim", obj))
rep <- sdreport(obj)

fit <- lmer(Y ~ x + (x|fgps), data=df, REML=FALSE)

##Random Effects from both models
TMB.rand <- as.matrix(rep$par.random)
TMB.rand.int <- TMB.rand[1:20,]
TMB.rand.slope <- TMB.rand[21:40,]
TMB.rand.col <- cbind(TMB.rand.int, TMB.rand.slope)

lmer.rand<-ranef(fit)

##plots
plot(ranef(fit)[[1]][,1], TMB.rand.int); abline(c(0,1))
plot(ranef(fit)[[1]][,2], TMB.rand.slope); abline(c(0,1))

