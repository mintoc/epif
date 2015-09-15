##---------------------------
## Various ways of fitting multinomial
## CM, BB: Tue Sep 15 2015
## roll own: multinomial with regressors in R
## Dirichlet multinomial
##---------------------------

## make some multinomial data
library(MASS)
library(nnet)

Y <- t(rmultinom(3, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))

multinom.fit <- multinom(Y ~ 1)

predict(multinom.fit, type = "prob")[1,]

## same as
colSums(Y)/sum(Y)

## fit own in R

X <- matrix(1, nrow = nrow(Y))
beta <- matrix(c(0, coef(multinom.fit)))
eta <- X %*% t(beta)
exp(eta) / (1 + sum(exp(eta[1, - 1])))

## write our own negative log-likelihood

multinom.nll <- function(theta, Y){
  beta <- matrix(c(0, theta)) ## note zero, a la Greene
  X <- matrix(1, nrow = nrow(Y))
  eta <- X %*% t(beta)
  P <- exp(eta) / (1 + sum(exp(eta[1, - 1])))
  ## note doing this because dmultinom works on a vector not matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    dmultinom(Y[z,], prob = P[z,], log = TRUE)
  }))
  return(nll)
}

## start at starting values c(0,0,0)
fit <- optim(par = c(0,0,0), fn = multinom.nll, Y = Y)

fit$par
coef(multinom.fit)
mglm.fit$coefficients

##
library(MGLM)
mglm.fit <- MGLMreg(Y ~ 1, dist = "MN")

##
library(mlogit)

Y.df0 <- data.frame(grp = rep(letters[1:4], each = 3), count = c(Y))

Y.df0 <- Y.df0[rep(1:12, times = Y.df0$count), ]

rownames(Y.df0) <- NULL

mlogit.dat <- mlogit.data(Y.df0, shape="wide", choice="grp")

mlogit.fit <- mlogit(grp ~ 0 | 1, data = mlogit.dat)

P <- mglm.fit$fitted

dmultinom(Y, prob = P, size = 20, log = TRUE)

log(factorial(20)/(prod(factorial(Y[1,]))) * prod(P[1,] ^ Y[1,]) * 
  factorial(20)/(prod(factorial(Y[2,]))) * prod(P[2,] ^ Y[2,]) *
  factorial(20)/(prod(factorial(Y[3,]))) * prod(P[3,] ^ Y[3,]))

## INCLUDE A REGRESSOR
Y <- t(rmultinom(20, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))

x1 <- rnorm(nrow(Y))
x2 <- rnorm(nrow(Y))
multinom.fit <- multinom(Y ~ x1 + x2)

multinom.nll <- function(theta, Y, X){
  npar <- (ncol(Y) - 1) *  ncol(X)
  if(length(theta) != npar){
    stop("Number of parameters supplied incorrect")
  }
  ##
  beta <- cbind(0, matrix(theta, ncol = ncol(Y) - 1)) ## note zero, a la Greene
  ## linear predictor
  eta <- X %*% beta
  ## probabilities - double check
  P <- exp(eta) / (1 + rowSums(exp(eta[, - 1])))
  ## log-likelihood
  ## note doing this because dmultinom works on a vector not matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    dmultinom(Y[z,], prob = P[z,], log = TRUE)
  }))
  return(nll)
}

X <- cbind(1, x1, x2)

npar <- (ncol(Y) - 1) *  ncol(X)
start.par <- rep(0, npar)

fit <- optim(par = start.par, fn = multinom.nll, Y = Y, X = X, control = list(trace = 1), method = "BFGS")

mglm.fit <- MGLMreg(Y ~ x1 + x2, dist = "MN")

## categorical covariate also
pet <- sample(c("terrapin","none","goldfish"), replace = TRUE, size = nrow(Y))

## how to grab a model matrix
multinom.fit <- multinom(Y ~ x1 + x2 + pet)
X <- model.matrix(multinom.fit)

npar <- (ncol(Y) - 1) *  ncol(X) 
start.par <- rep(0, npar)

fit <- optim(par = start.par, fn = multinom.nll, Y = Y, X = X, control = list(trace = 1), method = "BFGS")

##-----------------------
## DIRICHLET-MULTINOMIAL
##-----------------------
## alpha <- c(4, 3, 2,2)
## n <- c(1, 2, 6, 5)
## probability mass function for Dirichlet-multinomial

ddm <- function(n, alpha, log = FALSE){
  ## calculates pdf on log-scale
  ## to deal with large numbers
  ## http://www2.math.su.se/matstat/reports/seriec/2014/rep6/report.pdf
  N <- sum(n)
  lpdf <- lfactorial(N) - sum(lfactorial(n)) + lgamma(sum(alpha)) - lgamma(sum(alpha + n)) + sum(lgamma(alpha + n) - lgamma(alpha))
  pdf <- exp(lpdf)
  if(log){
    return(lpdf)
  }else{
    return(pdf)
  }
}

n <- c(10, 10, 20, 5)
alpha <- c(1, 2, 3, 1)

ddm(n =  n, alpha = alpha)

ddm2 <- function(n, alpha, log = FALSE){
  ## from wikipedia
  ## calculates pdf on log-scale
  ## to deal with large numbers 
  N <- sum(n)
  A <- sum(alpha)
  ##pdf <- gamma(A) / gamma(N + A) * prod(gamma(n + alpha) / gamma(alpha))
  lpdf <- lgamma(A) - lgamma(N + A) + sum(lgamma(n + alpha) - lgamma(alpha))
  pdf <- exp(lpdf)
  if(log){
    return(lpdf)
  }else{
    return(pdf)
  }
}

ddm2(n = n, alpha = alpha, log = TRUE)

## n1 <- n * 100
## ddm(n = n1, alpha = alpha)

dm.nll <- function(theta, Y){
  alpha <- exp(theta)
  ## note doing this because dmultinom works on a vector not matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    ddm2(n = Y[z,], alpha = alpha, log = TRUE)
  }))
  return(nll)
}

library(dirmult)
Y <- t(rmultinom(10, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))
##Y <- rdirm(size = 20, alpha = c(1,2,3,4), n = 10)

##
fit <- optim(par = rep(0, 4), fn = dm.nll, Y = Y, control = list(trace = 1), method = "BFGS")

exp(fit$par) / sum(exp(fit$par))

mglm.fit <- MGLMreg(Y ~ 1, dist = "DM")

mglm.fit$fitted[1, ]

##---------------------------------------
## DIRICHLET-MULTINOMIAL WITH REGRESSORS 
##---------------------------------------

dm.nll <- function(theta, Y, X){
  ## note using ddm not ddm2
  npar <- ncol(Y) *  ncol(X)
  if(length(theta) != npar){
    stop("Number of parameters supplied incorrect")
  }
  ##
  beta <- matrix(theta, ncol = ncol(Y))
  ## linear predictor
  eta <- X %*% beta
  ## alphas
  alpha <- exp(eta) ## see http://www2.math.su.se/matstat/reports/seriec/2014/rep6/report.pdf
  ## log-likelihood
  ## note doing this because dmultinom works on a vector not matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    ddm(n = Y[z,], alpha = alpha[z,], log = TRUE)
  }))
  return(nll)
}

Y <- t(rmultinom(20, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))

x1 <- rnorm(nrow(Y))
x2 <- rnorm(nrow(Y))
multinom.fit <- multinom(Y ~ x1 + x2)
X <- model.matrix(multinom.fit)

fit <- optim(par = rep(0, 12), fn = dm.nll, Y = Y, X = X, control = list(trace = 1, maxit = 1e4), method = "BFGS")

theta <- fit$par
beta <- matrix(theta, ncol = ncol(Y))
eta <- X %*% beta
alpha <- exp(eta) ## see http://www2.math.su.se/matstat/reports/seriec/2014/rep6/report.pdf
p.dm <- alpha / rowSums(alpha)

## in MGLM
mglm.fit <- MGLMreg(Y ~ x1 + x2, dist = "DM")

plot(mglm.fit$fitted, p.dm)
abline(c(0,1))

## NEXT STEP - FIT MULTINOMIAL IN TMB, NO RANDOM EFFECTS YET


##-----
## TMB 
##-----

## devtools::install_github("kaskr/adcomp", subdir = "TMB")
library(TMB)

## from https://github.com/kaskr/adcomp/wiki/Tutorial

n = 100
x = rnorm(n=n,mean=0,sd=1)          # Generate data

compile("tutorial.cpp")               # Compile the C++ file
dyn.load(dynlib("tutorial"))               # Dynamically link the C++ code

f = MakeADFun(data=list(x=x),parameters=list(mu=0,sigma=1))

f$fn(list(mu=1,sigma=1))          # Call TMB function value

-sum(dnorm(x,mean=1,sd=1,log=T))  # Verify that we get the same in R

f$gr(list(mu=1,sigma=1))  # First order derivatives (w.r.t mu and sigma)

f$he(list(mu=1,sigma=1))          # Second order derivatives, i.e. the Hessian matrix

fit = nlminb(f$par,f$fn,f$gr,lower=c(-10,0.0),upper=c(10.0,10.0))
print(fit)
