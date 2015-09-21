##---------------------------
## Various ways of fitting multinomial
## CM, BB: Tue Sep 15 2015
## roll own: multinomial with regressors in R
## Dirichlet multinomial
## Notes: using BFGS optimizer gets closer fits 
## to MGLM than Nelder-Mead
##---------------------------

library(MASS)
library(nnet)
library(MGLM)

## load the functions
source("multi_funs.R")

## MULTINOMIAL WITH REGRESSORS
## make some multinomial data
Y <- t(rmultinom(20, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))

## random regressors
x1 <- rnorm(nrow(Y))
x2 <- rnorm(nrow(Y))
## fit with nnet
multinom.fit <- multinom(Y ~ x1 + x2)

## fit with mglm
mglm.fit <- MGLMreg(Y ~ x1 + x2, dist = "MN")

## fit with own
X <- cbind(1, x1, x2)

npar <- (ncol(Y) - 1) *  ncol(X)
start.par <- rep(0, npar)

fit <- optim(par = start.par, fn = multinom.nll, Y = Y, X = X, control = list(trace = 1), method = "BFGS")

## multinomial fit
multinom.fit <- multinom(Y ~ x1 + x2)

## compare log-likelihoods
- fit$value
mglm.fit$logL
logLik(multinom.fit) ## different as multinomial coefficient ommitted here

## compare predictions
## own
beta <- cbind(0, matrix(fit$par, ncol = ncol(Y) - 1)) ## note zero, a la Greene
eta <- X %*% beta
P.own <- exp(eta) / (1 + rowSums(exp(eta[, - 1])))

## multinom
P.multinom <- predict(multinom.fit, type = "prob")

## mglm
P.mglm <- mglm.fit$fitted

P.own / P.multinom ## same
P.own / P.mglm ## close

## ILLUSTRATION WITH CATEGORICAL COVARIATE
pet <- sample(c("terrapin","none","goldfish"), replace = TRUE, size = nrow(Y))

## how to grab a model matrix
multinom.fit <- multinom(Y ~ x1 + x2 + pet)
X <- model.matrix(multinom.fit)

## use this in own fitting code
npar <- (ncol(Y) - 1) *  ncol(X) 
start.par <- rep(0, npar)

fit <- optim(par = start.par, fn = multinom.nll, Y = Y, X = X, control = list(trace = 1), method = "BFGS")

beta <- cbind(0, matrix(fit$par, ncol = ncol(Y) - 1)) ## note zero, a la Greene
eta <- X %*% beta
P.own <- exp(eta) / (1 + rowSums(exp(eta[, - 1])))

## multinom
P.multinom <- predict(multinom.fit, type = "prob")

P.own / P.multinom

plot(P.own, P.multinom)

##-----------------------
## DIRICHLET-MULTINOMIAL
##-----------------------
## alpha <- c(4, 3, 2,2)
## n <- c(1, 2, 6, 5)
## probability mass function for Dirichlet-multinomial
## NOTE: straightforward to include offset in own code

Y <- t(rmultinom(20, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))

x1 <- rnorm(nrow(Y))
x2 <- rnorm(nrow(Y))
X <- cbind(1, x1, x2)

fit <- optim(par = rep(0, 12), fn = dm.nll, Y = Y, X = X, control = list(trace = 1, maxit = 1e4), method = "BFGS")

## compare log-likelihoods
- fit$value
mglm.fit$logL

theta <- fit$par
beta <- matrix(theta, ncol = ncol(Y))
eta <- X %*% beta
alpha <- exp(eta) ## see http://www2.math.su.se/matstat/reports/seriec/2014/rep6/report.pdf
P.own <- alpha / rowSums(alpha)

## in MGLM
mglm.fit <- MGLMreg(Y ~ x1 + x2, dist = "DM")

P.mglm <- mglm.fit$fitted

P.own / P.mglm

plot(P.own, P.mglm)
abline(c(0,1))

##------
## ADMB
##------
## write the data out
datfile <- "../admb/multinomial.dat"
cat("# number of observations n \n", nrow(Y), "\n", file = datfile)
cat("# number of groups m \n", ncol(Y), "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", ncol(X), "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Y, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)


##-----
## TMB 
##-----
## NEXT STEP - FIT MULTINOMIAL IN TMB, NO RANDOM EFFECTS YET

Y <- t(rmultinom(20, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))
x1 <- rnorm(nrow(Y))
x2 <- rnorm(nrow(Y))
X <- cbind(1, x1, x2)

## devtools::install_github("kaskr/adcomp", subdir = "TMB")
library(TMB)

compile("multinomial.cpp")               # Compile the C++ file

dyn.load(dynlib("multinomial"))               # Dynamically link the C++ code
## dyn.unload(dynlib("multinomial"))

##f = MakeADFun(data=list(x=x), parameters=list(mu=0,sigma=1))

f <- MakeADFun(data=list(Y= Y, X = X), parameters = list(theta = rep(1, 9)))



f$fn(list(mu=1,sigma=1))          # Call TMB function value

-sum(dnorm(x,mean=1,sd=1,log=T))  # Verify that we get the same in R

f$gr(list(mu=1,sigma=1))  # First order derivatives (w.r.t mu and sigma)

f$he(list(mu=1,sigma=1))          # Second order derivatives, i.e. the Hessian matrix

fit = nlminb(f$par,f$fn,f$gr,lower=c(-10,0.0),upper=c(10.0,10.0))
print(fit)
