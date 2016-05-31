setwd("C:/Users/bburke/Desktop/GMIT/Gear Trials/R/TMB")
library(TMB)
library(Rcpp)
library(lme4)

tmb_rand<-"
#include<TMB.hpp>
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
  PARAMETER_VECTOR(a); // random
  PARAMETER_VECTOR(b); // random
  PARAMETER(logsigma_b); // optimize on the real line and then transform
  PARAMETER(logsigma_n);
  PARAMETER(logsigma_a);
//PROCEDURE:
  Type sigma_b = exp(logsigma_b);
  Type sigma_n = exp(logsigma_n);
  Type sigma_a = exp(logsigma_a);
//
Type nll = 0.0;
for(int j = 0; j < m; j++){
  nll -= dnorm(a(j), Type(0), sigma_a, true);
  nll -= dnorm(b(j), Type(0), sigma_b, true);
}
//
for(int i = 0; i < n; i++){
  int indexi = gp[i];
    nll -= dnorm(Y[i], (alpha + a[indexi]) + (beta + b[indexi])*X[i], sigma_n, true);
}
return nll;
}
"

write(tmb_rand, file = "rand_model.cpp")

compile("rand_model.cpp")

dyn.unload(dynlib("rand_model")) ## before loading to be safe
dyn.load(dynlib("rand_model"))

rm(obj)

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

library(ggplot2)
ggplot(df, aes(x=x, y=y)) + geom_point() + facet_wrap(~gps)

Y <- matrix(df$y, ncol = 1)

## note random statement

obj <- MakeADFun(
  data = list(Y = df$y, X = df$x , n = nrow(Y), m = m, gp = df$gps),
  parameters = list(alpha = 0, beta = 0, a = rep(0,m), b = rep(0,m), logsigma_b = 0, logsigma_n = 0, logsigma_a = 0),
  random = c("a","b"),
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

