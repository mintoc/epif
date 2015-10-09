library(MASS)
library(nnet)

## example where true p = 0.25 for each case but they are subsampled differently
## observed counts
Y <- cbind(10,20,8,40)
## sub-sampling ratios
q <- cbind(0.125,0.25,0.1,0.5)

## naive fit without subsampling
fit0 <- multinom(Y~1)
predict(fit0, type = "prob")
## same as ratio of raw counts
Y/sum(Y)

## make a data frame for the fit with the offset
offset.mat <- log(q/q[1])
df <- data.frame(Y,offset.mat)
names(df) <- c(paste("y",1:4, sep = ""), paste("off",1:4, sep = ""))

fit1 <- multinom(cbind(y1,y2,y3,y4) ~ 1 + offset(cbind(off1, off2, off3, off4)), data = df)

## predict for equal sampling
predict(fit1, type = "prob", newdata = data.frame(off1 = 0, off2 = 0, off3 = 0, off4 = 0))

## OWN CODE

## MULTINOMIAL NEGATIVE LOG-LIKELIHOOD
multinom.offset.nll <- function(theta, Y, X, Off){
  npar <- (ncol(Y) - 1) *  ncol(X)
  if(length(theta) != npar){
    stop("Number of parameters supplied incorrect")
  }
  ##
  beta <- cbind(0, matrix(theta, ncol = ncol(Y) - 1)) ## note zero, a la Greene
  ## linear predictor
  eta <- X %*% beta + Off
  ## probabilities - double check with Greene
  ##P <- exp(eta) / (1 + rowSums(exp(eta[, - 1])))
  ##P <- exp(eta) / (1 + rowSums(exp(eta)))
  P <- exp(eta) / rowSums(exp(eta))
  ## log-likelihood
  ## note doing this because dmultinom works on a vector not matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    dmultinom(Y[z,], prob = P[z,], log = TRUE)
  }))
  return(nll)
}

## fit with own
X <- model.matrix(fit1)
npar <- (ncol(Y) - 1) *  ncol(X)
start.par <- rep(0, npar)
Off <- offset.mat

fit <- optim(par = start.par, fn = multinom.offset.nll, Y = Y, X = X, Off = Off, control = list(trace = 1), method = "BFGS")

## get predicted probability
beta <- cbind(0, matrix(fit$par, ncol = ncol(Y) - 1)) ## note zero, a la Greene
eta <- X %*% beta ## offsets all zero here
P <- exp(eta) / rowSums(exp(eta))

