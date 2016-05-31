
library(TMB)
library(Rcpp)

m <- 20
n <- 10
mu <- 2
sigma.m <- 3 ## random effects SD
sigma.n <- 1 ## residual SD
obs <- m * n

b <- rnorm(m, mean = mu, sd = sigma.m)

y.mat <- matrix(NA, nrow = n, ncol = m)

for(i in 1:m){
  y.mat[, i] <- rnorm(n, mean = b[i], sd = sigma.n)
}

df <- data.frame(gps = rep(0:(m-1),each = n), y = c(y.mat))

tmb_rand<-"
#include<TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()(){

//DATA:
DATA_VECTOR(Y);
DATA_IVECTOR(gp);
DATA_SCALAR(m);//groups
DATA_SCALAR(n);//individuals

//PARAMETER:
PARAMETER(mu);
PARAMETER_VECTOR(b); // random 
PARAMETER(logsigma_m); // optimize on the real line and then transform
PARAMETER(logsigma_n);

//PROCEDURE:
Type sigma_m = exp(logsigma_m);
Type sigma_n = exp(logsigma_n);

Type nll = 0.0;

for(int i = 0; i < m; i++){
  nll -= dnorm(b(i), Type(0), sigma_m, true);
}

for(int i = 0; i < n; i++){
  int indexi = gp[i];
  nll -= dnorm(Y[i], mu + b[indexi], sigma_n, true);
}
return nll;
}
"

write(tmb_rand, file = "rand_model.cpp")

compile("rand_model.cpp")

dyn.unload(dynlib("rand_model")) ## before loading to be safe
dyn.load(dynlib("rand_model"))

## strangely need to remove obj here if changing the data
rm(obj)

Y <- matrix(df$y, ncol = 1)

## note random statement

obj <- MakeADFun(
         data = list(Y = df$y, n = nrow(Y), m = m, gp = df$gps),
         parameters = list(mu = 0, b = rep(0,m), logsigma_m = 0, logsigma_n = 0),
         random = c("b"),
         DLL = "rand_model")

(opt<-do.call("optim", obj))
rep <- sdreport(obj)

bhat <- summary(rep, "random")[, "Estimate"]
pred <- rep$par.fixed["mu"] + bhat
pred.df <- data.frame(gps = 0:(m - 1), y = pred)

## 
library(lme4)
lme.fit <- lmer(y ~ 1 + (1|gps), data = df, REML = FALSE)

plot(bhat, ranef(lme.fit)[[1]][, 1]); abline(c(0,1))

library(ggplot2)
theme_set(theme_bw())

ggplot(df, aes(x = gps, y = y)) + geom_point(colour = "slategrey") +
  geom_point(data = pred.df, colour = "red", aes(size = 2)) + theme(legend.position = "none")
