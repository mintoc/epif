##PDF
ddm <- function(n, alpha, log = FALSE){
  ## calculates pdf on log-scale
  ## to deal with large numbers
  ## http://www2.math.su.se/matstat/reports/seriec/2014/rep6/report.pdf
  N <- sum(n)
  n <- matrix(n, nrow = 1)
  alpha <- matrix(alpha, nrow = 1)
  ## rowSums version
  lpdf <- lgamma(N + 1) - rowSums(lgamma(n+1)) +
          lgamma(rowSums(alpha)) - lgamma(rowSums(alpha + n)) + 
          rowSums(lgamma(alpha + n) - lgamma(alpha))
  pdf <- exp(lpdf)
  if(log){
    return(lpdf)
  }else{
    return(pdf)
  }
}


###Log Likelihood Function - Offset Ommitted
dm.nll <- function(theta, Y, X){
    print(theta)
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


###Log Likelihood Function - Offset Included
dm.nll.o <- function(theta, Y, X, off){
  print(theta)
  ## note using ddm not ddm2
  npar <- ncol(Y) *  ncol(X)
  if(length(theta) != npar){
    stop("Number of parameters supplied incorrect")
  }
  ##
  beta <- matrix(theta, ncol = ncol(Y))
  ## linear predictor
  eta <- X %*% beta + log(off)
  ## alphas
  alpha <- exp(eta) ## see http://www2.math.su.se/matstat/reports/seriec/2014/rep6/report.pdf
  ## log-likelihood
  ## note doing this because ddm works on a vector not a matrix
  nll <- - sum(sapply(1:nrow(Y), function(z){
    ddm(n = Y[z,], alpha = alpha[z,], log = TRUE)
  }))
  return(nll)
}


###Gradient Function
objfun.grad <- function(alpha, x, y, d, p){
  alpha <- matrix(alpha, p, d)
  Beta <- exp(x%*%alpha)
  m <- rowSums(y)
  tmpvector <- digamma(rowSums(Beta)+m)-digamma(rowSums(Beta))
  tmpvector[is.nan(tmpvector)] <- 0
  tmpmatrix <- digamma(Beta + y) - digamma(Beta)
  tmpvector2 <- trigamma(rowSums(Beta)) - trigamma(m+rowSums(Beta))
  tmpmatrix2 <- trigamma(Beta) - trigamma(Beta+y)
  tmpmatrix2 <- Beta*tmpmatrix - Beta^2*tmpmatrix2
  dalpha <- Beta*(tmpmatrix - tmpvector)
  
  expr <- paste("rbind(", paste(rep("dalpha", p), collapse = ","), 
                ")", sep = "")
  dalpha <- eval(parse(text=expr))
  dalpha <- matrix(c(dalpha), nrow(x), ncol=p*d)
  expr2 <- paste("cbind(", paste(rep("x", d), collapse = ","), 
                 ")", sep = "")
  x <- eval(parse(text = expr2))
  dl <- colSums(dalpha*x)
  return(-dl)
}
