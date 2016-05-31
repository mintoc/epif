##------------------------------------------
## Simulated example of multinomial mixed effects 
## CM: Thu Dec 3 2015
## Simulates data similar to catch-comparison trial for Nephrops
##------------------------------------------

set.seed(101)

## number of cod-ends
m <- 4

## number of hauls (clusters)
ngp <- 15

## carapace length (continuous covariate)
cl <- seq(15, 45, length = 30)

dat <- expand.grid(cl = cl, haul = paste("H", 1:ngp, sep = ""))

## net position configuration
dat$netpos <- ifelse(dat$haul %in% c("H1", "H2", "H3"), "np1",
                     ifelse(dat$haul %in% c("H4", "H5", "H6"), "np2",
                            ifelse(dat$haul %in% c("H7", "H8", "H9"), "np3", "np4")))

## multinomial logit model matrix
dat$y <- rnorm(nrow(dat))
ff <- y ~ cl + netpos
str(mm <- model.frame(ff, dat))
X <- model.matrix(ff, mm)
n <- nrow(X)

## conditional logit model matrix, e.g., total cod-end weight
Xcond0 <- matrix(rnorm(ngp * m, mean = 350, sd = 50), ncol = m)
Xcond <- Xcond0[rep(1:ngp, each = length(cl)), ]

## fixed effects multinomial logit coefficients
beta <- cbind(0,
              rbind(runif(m-1, min = -2, max = 0), ## intercept
                    runif(m-1, min = 0, max = 0.04), ## cl effect
                    matrix(runif((m-1)^2, min = -1, max = 1), ncol = m-1) ## net position effects
                    )
              )

## fixed effects conditional logit coefficient
betacond <- rnorm(1, sd = 1e-3)

## random effects
library(MASS)
## vcov
##a <- rnorm(6, sd = 0.3) ## working
a <- runif(6, min = 0.1, max = 0.4)
L <- matrix(0, ncol = m - 1, nrow = m - 1)

ii <- 1
for(i in 1:(m-1)){
  for(j in 1:i){
    L[i, j] <- a[ii]
    ii <- ii + 1
  }
}
Sigma <- L %*% t(L)
u.mat0 <- mvrnorm(ngp, mu = rep(0, m - 1), Sigma = Sigma)

plot(as.data.frame(u.mat0)) ## check not too strongly correlated, if so, change seed

u.mat <- cbind(0, u.mat0[rep(1:ngp, each = length(cl)), ])

gps <- as.numeric(dat$haul)

## linear predictor
eta <- X %*% beta + Xcond * betacond + u.mat

p <- exp(eta) / rowSums(exp(eta))

plot.dat <- data.frame(p = c(p), cl = rep(dat$cl, m), haul = rep(dat$haul, m), mesh = rep(paste("mesh", 1:m, sep = ""), each = nrow(X)))

library(ggplot2)
ggplot(plot.dat, aes(x = cl, y = p)) + geom_line() + facet_grid(haul ~ mesh)

## response
Yall <- t(apply(p, 1, FUN = function(z){rmultinom(1, size = 100, prob = z)}))

## subsampling 
subs.mat0 <- matrix(runif(ngp * m), ncol = m)
##subs.mat0 <- matrix(rep(1, ngp * m), ncol = m)
subs.mat <- subs.mat0[rep(1:ngp, each = length(cl)), ]

Yobs <- round(Yall * subs.mat)

Offset <- log(subs.mat / subs.mat[,1])

## prediction matrices
Xpred <- cbind(X[, 1:2], 1/3, 1/3, 1/3)
Xcondpred <- Xcond * 0

npred <- nrow(Xpred)

## send to admbre
## write the data out
datfile <- "../admbre/multinomialme.dat"
cat("# number of observations n \n", n, "\n", file = datfile)
cat("# number of categories m \n", m, "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", ncol(X), "\n", file = datfile, append = TRUE)
cat("# dimension of conditional variables q \n", ncol(Xcond), "\n", file = datfile, append = TRUE)
cat("# number of groups ngp \n", ngp, "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Yobs, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# conditional matrix Xcond \n", file = datfile, append = TRUE)
write.table(Xcond, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# groups \n", gps, "\n", file = datfile, append = TRUE)
cat("# Offset (log(q[i]/q[1])) \n", file = datfile, append = TRUE)
write.table(Offset, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
## predictions
cat("# number of prediction rows npred \n", npred, "\n", file = datfile, append = TRUE)
cat("# prediction matrix Xpred \n", file = datfile, append = TRUE)
write.table(Xpred, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# conditional prediction matrix Xcondpred \n", file = datfile, append = TRUE)
write.table(Xcondpred, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## starting values pinfile
beta0.start <- matrix(0, nrow = nrow(beta), ncol = ncol(beta) - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)
nchol <- (m-1)*(m)/2
L.start <- diag(m-1)
L.start[upper.tri(L.start)] <- NA
a.start <- as.numeric(na.omit(c(t(L.start))))
##a.start <- rep(0.1, 3)

pinfile <- "../admbre/multinomialme.pin"
cat("# beta0 \n", file = pinfile)
write.table(beta0.start, file = pinfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# betacond \n", 0, "\n", file = pinfile, append = TRUE)
##set.seed(10)
cat("# a \n", a.start, "\n", file = pinfile, append = TRUE)
##cat("# a \n", rep(0.1, m-1), "\n", file = pinfile, append = TRUE)
cat("# u0 \n", file = pinfile, append = TRUE)
write.table(u0.start, file = pinfile, append = TRUE, col.names = FALSE, row.names = FALSE)



## RUN THE MODEL IN ADMB-RE - CODE BELOW WON'T WORK OTHERWISE

## read in and plot the results
coef.admbre <- read.table("../admbre/multinomialme.std", header = TRUE)

uhat <- matrix(subset(coef.admbre, name == "u0")$value, ncol = m-1, byrow = TRUE)

plot(u.mat0, uhat)
abline(c(0, 1)) ## not very well recovered in this example

## Covariation matrix of the random effects
a.hat <- subset(coef.admbre, name == "a")$value

ii <- 1

L <- matrix(0, 3, 3)
for(i in 1:3){
  for(j in 1:i){
    L[i,j] <- a.hat[ii]
    ii <- ii + 1
  }
}

sigma.hat <- L %*% t(L)

plot(Sigma, sigma.hat); abline(c(0,1))

##
beta.hat <- cbind(0, matrix(subset(coef.admbre, name == "beta0")$value, nrow = ncol(X), byrow = TRUE))

plot(beta.hat[, -1], beta[, -1])
abline(c(0, 1))
