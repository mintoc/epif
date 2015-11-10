##-------------------
## Quad-rig analyses
## CM, BB: Tue Sep 22 2015
## Note: bigger models run on the cluster
## go through each combination outputting to admb and
## running on the cluster
##-------------------

library(nnet)

## load the data
load("our_lass_data_objects.RData")

## load additional functions
##source("funs.R")

## using poly otherwise very strong correlations induced
mnom.fit <- multinom(neph.count.mat ~ 
                     Carapace.length + 
                     netconfig +
                     offset(log(cbind(m70mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m80mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m90mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m100mm_SUBSRATIO/m70mm_SUBSRATIO))),
                     data = our.lass.neph.cast)

## pred.P <- predict(mnom.fit, type = "prob")
## Y.raised <- neph.count.mat / subsratio.mat
## pred.Y <- apply(pred.P, 2, "*", rowSums(Y.raised))
## plot(pred.Y, Y.raised, col = "grey"); abline(c(0,1))

## mnom.resid <- residuals(mnom.fit)

## plot.multinom.resid("Carapace.length", ylim = c(-.5, .5))
## plot.multinom.resid("netconfig", ylim = c(-.5, .5))

## plot.multinom.resid("m70mm_Total.catch", ylim = c(-.5, .5))
## plot.multinom.resid("m80mm_Total.catch", ylim = c(-.5, .5))
## plot.multinom.resid("m90mm_Total.catch", ylim = c(-.5, .5))
## plot.multinom.resid("m100mm_Total.catch", ylim = c(-.5, .5))

## library(fields)
## par.corr <- cov2cor(vcov(mnom.fit))
## image.plot(par.corr)
## hist(par.corr[lower.tri(par.corr)], breaks = 100, xlim = c(-1, 1), xaxs = "i")

## ADMB-RE
Y <- neph.count.mat
X <- model.matrix(mnom.fit)
gps <- as.numeric(our.lass.neph.cast$fHAUL)
ngp <- length(unique(gps))

## Conditional variable - catch per net
Xcond <- our.lass.neph.cast[, c("m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")]

## predictions matrix
min.cl <- min(X[,"Carapace.length"])
max.cl <- max(X[,"Carapace.length"])
cl.vec <- seq(min.cl, max.cl, length = 100)

##head(model.matrix(mnom.fit))
##Xpred <- model.matrix(mnom.fit)[1:2,]
Xpred <- cbind(1, cl.vec, 1/3, 1/3, 1/3)

##mean.bulk <- mean(unlist(bulk.cast[, -1])) ## 372.599
mean.bulk <- 372.599
Xcondpred <- Xcond[1,]
Xcondpred[] <- mean.bulk
Xcondpred <- Xcondpred[rep(1, nrow(Xpred)),]
npred <- nrow(Xcondpred)

## NOW SEND THESE OVER TO ADMB
n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)
q <- ncol(Xcond)
nchol <- (m-1)*(m)/2

beta0.start <- matrix(0, nrow = p, ncol = m - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)

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
## predictions
cat("# number of categories npred \n", npred, "\n", file = datfile, append = TRUE)
cat("# prediction matrix Xpred \n", file = datfile, append = TRUE)
write.table(Xpred, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# conditional prediction matrix Xcondpred \n", file = datfile, append = TRUE)
write.table(Xcondpred, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## pinfile
datfile <- "../admb/conlogit_re/conditional.pin"
cat("# beta0 \n", file = datfile)
write.table(beta0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# betacond \n", 0, "\n", file = datfile, append = TRUE)
##set.seed(10)
cat("# a \n", rep(1, nchol), "\n", file = datfile, append = TRUE)
##cat("# a \n", rep(0.1, m-1), "\n", file = datfile, append = TRUE)
cat("# u0 \n", file = datfile, append = TRUE)
write.table(u0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## SAVE OBJECTS
save(list = c("Y", "X", "Xpred"), file = "model_objects.RData")

## run e vous


##-----------------------------
## Dirichlet-Multinomial model
##-----------------------------


