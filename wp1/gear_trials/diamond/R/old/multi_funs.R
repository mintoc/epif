## Functions for fitting multi-component count data

## MULTINOMIAL NEGATIVE LOG-LIKELIHOOD
multinom.nll <- function(theta, Y, X){
  npar <- (ncol(Y) - 1) *  ncol(X)
  if(length(theta) != npar){
    stop("Number of parameters supplied incorrect")
  }
  ##
  beta <- cbind(0, matrix(theta, ncol = ncol(Y) - 1)) ## note zero, a la Greene
  ## linear predictor
  eta <- X %*% beta
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

## DIRICHLET-MULTINOMIAL PDF
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

## DIRICHLET-MULTINOMIAL NEGATIVE LOG-LIKELIHOOD
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
  ## note doing this because ddm works on a vector not a matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    ddm(n = Y[z,], alpha = alpha[z,], log = TRUE)
  }))
  return(nll)
}
