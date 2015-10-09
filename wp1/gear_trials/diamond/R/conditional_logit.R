##--------------------------
## Conditional logit trials 
## CM:
##
##--------------------------

library(mlogit)

data("Fishing", package = "mlogit")

Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")

## a pure "conditional" model
summary(mlogit(mode ~ price | 0, data = Fish))

## a pure "multinomial model"
fit <- mlogit(mode ~ 0 | 1, data = Fish)


## a pure "multinomial model"
summary(mlogit(mode ~ 0 | income, data = Fish))


## simple mixed example to work with
fit <- mlogit(mode ~ price | 1, data = Fish)


logLik(fit)

library("nnet")
summary(multinom(mode ~ 1, data = Fishing))


## a "mixed" model
m <- mlogit(mode ~ price+ catch | income, data = Fish)
summary(m)

## In ADMB
Fishing$mode
## ADMB-RE
beach boat charter pier 
Y <- with(Fishing,
          cbind(
            ifelse(mode == "beach", 1, 0),
            ifelse(mode == "boat", 1, 0),
            ifelse(mode == "charter", 1, 0),
            ifelse(mode == "pier", 1, 0)
            )
          )
n <- nrow(Y)
X <- matrix(1, nrow = n)
## Xcond formulated as difference from baseline
Xcond0 <- Fishing[, c("price.beach", "price.boat", "price.charter", "price.pier")]
##Xcond <- Xcond0[,2:4] - Xcond0[,1]
Xcond <- Xcond0[,1:4]

gps <- sample(1:5, n, replace = TRUE)
ngp <- length(unique(gps))

n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)
q <- ncol(Xcond)
nchol <- (m-1) * (m) / 2

beta0.start <- matrix(0, nrow = p, ncol = m - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)

offset.mat <- matrix(0, ncol = m, nrow = n)

## write the data out
datfile <- "../admb/conlogit_re/conditional.dat"
cat("# number of observations n \n", n, "\n", file = datfile)
cat("# number of categories m \n", m, "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", p, "\n", file = datfile, append = TRUE)
cat("# dimension of conditional variables q \n", q, "\n", file = datfile, append = TRUE)
cat("# number of groups ngp \n", ngp, "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Y, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# conditional matrix Xcond \n", file = datfile, append = TRUE)
write.table(Xcond, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# groups \n", gps, "\n", file = datfile, append = TRUE)
cat("# Offset (log(q[i]/q[1])) \n", file = datfile, append = TRUE)
write.table(offset.mat, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## pinfile
datfile <- "../admb/conlogit_re/conditional.pin"
cat("# beta0 \n", file = datfile)
write.table(beta0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# betacond \n", 0, "\n", file = datfile, append = TRUE)
##set.seed(10)
##cat("# a \n", rep(0.1, nchol), "\n", file = datfile, append = TRUE)
cat("# a \n", rep(0.1, m-1), "\n", file = datfile, append = TRUE)
cat("# u0 \n", file = datfile, append = TRUE)
write.table(u0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
