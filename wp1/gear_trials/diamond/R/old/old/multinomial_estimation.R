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
n <- 40
Y <- t(rmultinom(n, size = 20, prob = c(0.1, 0.3, 0.1, 0.5)))

## random regressors
x1 <- rnorm(nrow(Y))
x2 <- rnorm(nrow(Y))

## groups for random effects
gps <- rep(seq(1,n/4), each = 4)

## fit with nnet, no random effects
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
library(MASS)
## SIMULATE WITH EFFECTS
m <- 4
a <- rnorm((m-1)*(m)/2, sd = 1)

L <- matrix(0, nrow = m - 1, ncol = m - 1)

ii <- 1

for(i in 1:(m-1)){
  for(j in 1:i){
    L[i,j] <- a[ii]
    ii <- ii + 1
  }
}

Sigma <- L%*%t(L)

ngp <- 20
u0 <- mvrnorm(ngp, mu = rep(0, m-1), Sigma)
u <- cbind(0, u0)
nobspergp <- 10
n <- ngp * nobspergp

X <- cbind(1, rnorm(n))

p <- ncol(X)

beta0 <- matrix(rnorm((m-1) * p), nrow = p)

beta <- cbind(0, beta0)

gps <- rep(1:ngp, each = nobspergp)
eta <- X %*% beta + u[gps,]
P <- exp(eta) / (rowSums(exp(eta)))

Y <- t(apply(P, 1, FUN = function(z){
  rmultinom(1, size = 20, prob = z)
}))

## write the data out
datfile <- "../admb/multinomial_re/multinomial.dat"
cat("# number of observations n \n", nrow(Y), "\n", file = datfile)
cat("# number of categories m \n", ncol(Y), "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", ncol(X), "\n", file = datfile, append = TRUE)
cat("# number of groups ngp \n", length(unique(gps)), "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Y, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# groups \n", gps, "\n", file = datfile, append = TRUE)

## pinfile
datfile <- "../admb/multinomial_re/multinomial.pin"
cat("# beta0 \n", file = datfile)
write.table(beta0 * 0, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# a \n", rep(0.5, length(a)), "\n", file = datfile, append = TRUE)
cat("# u0 \n", file = datfile, append = TRUE)
write.table(u0 * 0, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

admb.res <- read.table("../admb/multinomial_re/multinomial.std", header = TRUE)

uhat.admb <- matrix(admb.res$value[admb.res$name == "u0"], ncol = m - 1, byrow = TRUE)

plot(u0, uhat.admb)
abline(c(0,1))

## BINARY EXAMPLE fit
library(lme4)
(gm1 <- glmer(Y ~ -1 + X + (1 | factor(gps)),
              family = binomial))

test <- multinom(Y ~ -1 + X)

logLik(gm1)

plot(ranef(gm1)[[1]][,1], uhat)
abline(c(0,1))


##------
## JAGS 
##------
## see: https://github.com/johnmyleswhite/JAGSExamples/blob/master/scripts/multinomial/multinomial.R
# Basic inference.
library(R2jags)

n <- nrow(Y)
k <- ncol(Y)
p <- ncol(X)
npar <- p * (m-1)
R <- diag(k-1)

# data
jags.data <- list("Y", "k", "n", "X", "p", "gps", "ngp", "R")
jags.params <- c("beta","sigma","u")

jags.inits <- function(){
  list("beta" = matrix(rnorm(npar), nrow = 2),
       "tau" = rWishart(1, df = k, Sigma = diag(k-1))[,,1],
       "u" = matrix(0,nrow = ngp, ncol = k-1)
       )
}

jagsfit <- jags(data = jags.data,
                inits = jags.inits,
                jags.params,
                n.iter = 1e4,
                model.file = "../jags/multinomial.bug")

plot(jagsfit)

print(jagsfit)

uhat.jags <- jagsfit$BUGSoutput$mean$u

## PLOT THEM ALL
true.df <- data.frame(method = "true",
                      u = c(u0),
                      k = rep(paste("u", 1:(m-1), sep = ""), each = ngp))

true.df$index <- seq(1, nrow(true.df))

admb.df <- data.frame(method = "admb",
                      u = c(uhat.admb),
                      k = rep(paste("u", 1:(m-1), sep = ""), each = ngp))
admb.df$index <- seq(1, nrow(admb.df))

jags.df <- data.frame(method = "jags",
                      u = c(uhat.jags),
                      k = rep(paste("u", 1:(m-1), sep = ""), each = ngp))
jags.df$index <- seq(1, nrow(jags.df))

all.df <- rbind(true.df, admb.df, jags.df)

library(ggplot2)

theme_set(theme_bw(base_size = 14))
ggplot(all.df, aes(x = index, y = u)) + geom_point(aes(colour = method, pch = method), size = 3) + facet_wrap(~k, scales = "free")


plot(u0, ujags)
abline(c(0,1))


