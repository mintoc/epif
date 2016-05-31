
## baseline counts:
N <- c(50, 25, 25)
prop.table(N)

## low weight
low.change <- c(0.7, 0.2, 0.4)
N.low <- N * low.change
p.low <- prop.table(N.low)

## high weight
high.change <- c(1.2, 1.5, 1.4)
N.high <- N * high.change
p.high <- prop.table(N.high)

get.p <- function(x){exp(x) / rowSums(exp(x))}

beta <- matrix(c(0, log(c(.25, .25) / .5)), nrow = 1)
eta <- matrix(1) %*% beta

get.p(eta)

cond.func <- function(theta, p){
    pred.p <- get.p(matrix(1) %*% beta + matrix(1, ncol = 3) %*% diag(theta))
    return(sum((p - pred.p)^2))
}

opt.low <- optim(fn = cond.func, par = c(1,1,1), p = prop.table(N.low), method = "BFGS", hessian = TRUE)
opt.high <- optim(fn = cond.func, par = c(1,1,1), p = prop.table(N.high), method = "BFGS")

(pred.low <- get.p(matrix(1) %*% beta + matrix(1, ncol = 3) %*% diag(opt.low$par)))
(pred.high <- get.p(matrix(1) %*% beta + matrix(1, ncol = 3) %*% diag(opt.high$par)))

## with zero in front
cond.func <- function(theta, p){
    pred.p <- get.p(matrix(1) %*% beta + matrix(1, ncol = 3) %*% diag(c(0, theta)))
    return(sum((p - pred.p)^2))
}

opt.low <- optim(fn = cond.func, par = c(1,1), p = prop.table(N.low), method = "BFGS", hessian = TRUE)
opt.high <- optim(fn = cond.func, par = c(1,1), p = prop.table(N.high), method = "BFGS", hessian = TRUE)

(pred.low <- get.p(matrix(1) %*% beta + matrix(1, ncol = 3) %*% diag(c(0, opt.low$par))))
(pred.high <- get.p(matrix(1) %*% beta + matrix(1, ncol = 3) %*% diag(c(0, opt.high$par))))

## generate some continuous data and see can you estimate?

## generate the counts
N0 <- c(50, 25, 25) * 1e3

weights <- matrix(0:100, nc = 1)

weight.rate <- matrix(c(-0.1, -0.05, -0.2), nr = 1)

N <- t(apply(exp(weights %*% weight.rate), 1, FUN = function(x){x * N0}))
n <- nrow(N)
p <- prop.table(N, 1)

matplot(p, type = "l", ylim = c(0,1), lty = 1)

fn <- function(theta){
    eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1)
    p <- get.p(eta)
    ll <- sapply(1:n, function(x){dmultinom(round(N[x,]), prob = p[x,], size = sum(round(N[x,])), log = TRUE)})
    return(-sum(ll))
}

matplot(p, type = "l", ylim = c(0,1), lty = 1)
matplot(prop.table(N, 1), type = "l", ylim = c(0,1), lty = 1)

opt <- optim(fn = fn, par = c(0,0))
theta <- opt$par
eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1)
p <- get.p(eta)

matlines(p)
## looks like it's not fitting well but the numbers out right are tiny

## now try to estimate effects of weight
fn <- function(theta){
    eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1) + weights[, rep(1,3)] %*% diag(theta[3:5])
    p <- get.p(eta)
    ll <- sapply(1:n, function(x){dmultinom(round(N[x,]), prob = p[x,], size = sum(round(N[x,])), log = TRUE)})
    return(-sum(ll))
}

matplot(prop.table(N, 1), type = "l", ylim = c(0,1), lty = 1)

opt <- optim(fn = fn, par = rep(0,5), hessian = TRUE)
theta <- opt$par
eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1) + weights[, rep(1,3)] %*% diag(theta[3:5])
pred.p <- get.p(eta)

matlines(pred.p, lwd = 2, lty = 2)

## hessian won't invert, can we estimate another way?
## have different weights in the nets

weights <- cbind(
             matrix(rnorm(100, mean = 50, sd = 10), nc = 1),
             matrix(rnorm(100, mean = 5, sd = 2), nc = 1),
             matrix(rnorm(100, mean = 25, sd = 5), nc = 1))

weight.rate <- diag(c(-0.1, -0.5, -0.08))

N <- t(apply(exp(weights %*% weight.rate), 1, FUN = function(x){x * N0}))
n <- nrow(N)
p <- prop.table(N, 1)

## now try to estimate effects of weight
fn <- function(theta){
    eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1) + weights %*% diag(c(0, theta[3:4]))
    p <- get.p(eta)
    ll <- sapply(1:n, function(x){dmultinom(round(N[x,]), prob = p[x,], size = sum(round(N[x,])), log = TRUE)})
    return(-sum(ll))
}

matplot(prop.table(N, 1), type = "l", ylim = c(0,1), lty = 1)

opt <- optim(fn = fn, par = rep(0, 4), hessian = TRUE, control = list(maxit = 1e4))
theta <- opt$par
eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1) + weights %*% diag(c(0, theta[3:4]))
pred.p <- get.p(eta)

matlines(pred.p, lwd = 2, lty = 2)

fn <- function(theta){
  eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1) + weights %*% diag(c(theta[3:5]))
  p <- get.p(eta)
  ll <- sapply(1:n, function(x){dmultinom(round(N[x,]), prob = p[x,], size = sum(round(N[x,])), log = TRUE)})
  return(-sum(ll))
}

matplot(prop.table(N, 1), type = "l", ylim = c(0,1), lty = 1)

opt <- optim(fn = fn, par = rep(0, 5), hessian = TRUE, control = list(maxit = 1e4))
theta <- opt$par
eta <- matrix(1, nr = n) %*% matrix(c(0, theta[1:2]), nr = 1) + weights %*% diag(c(theta[3:5]))
pred.p <- get.p(eta)

matlines(pred.p, lwd = 2, lty = 2)

## fits much better and hessian is invertible
