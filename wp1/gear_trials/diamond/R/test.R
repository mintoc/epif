

p <- colSums(Y)/sum(Y)

opt.alpha <- function(theta){
  alpha <- exp(theta)
  sum((alpha / sum(alpha) - p)^2)
}

fit <- optim(par = rep(1, 4), fn = opt.alpha)

alpha <- exp(fit$par)

p <- colSums(Y)/sum(Y)




