##-------------------
## Quad-rig analyses
## CM, BB: Tue Sep 22 2015
## Note: bigger models run on the cluster
## Run in TMB
## Also includes ADMB fits for comparison 
##-------------------
library(nnet)
library(TMB)

##------
## DATA
##------

load("our_lass_data_objects.RData")

## load additional functions
## source("funs.R")

scale.fun <- function(x){(x - mean(x)) / sd(x)}

## using poly otherwise very strong correlations induced
mnom.fit <- multinom(neph.count.mat ~
                     scale.fun(Carapace.length) + ##scale.fun(Carapace.length^2) +
                     ##netconfig +
                     offset(log(cbind(m70mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m80mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m90mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m100mm_SUBSRATIO/m70mm_SUBSRATIO))),
                     data = our.lass.neph.cast)

Y <- neph.count.mat
X <- model.matrix(mnom.fit)
n <- nrow(Y)
m <- ncol(Y)

gps <- as.numeric(our.lass.neph.cast$fHAUL)
ngp <- length(unique(gps))

##-----
## TMB 
##-----

library(TMB)

## model spec
tmb_multinomial_re<- "
  #include <TMB.hpp>
  using namespace density;
  template <class Type>
  Type objective_function<Type>::operator () (){
    //------
    // DATA
    //------
    DATA_MATRIX(Y);
    DATA_MATRIX(X);
    DATA_ARRAY(Xcond);
    DATA_MATRIX(Offset);
    DATA_IVECTOR(gp);
    DATA_INTEGER(ngp);
    DATA_MATRIX(W);
    DATA_MATRIX(WxCL);
    // PREDICTION
    DATA_MATRIX(Xpred);
    DATA_ARRAY(Xcondpred);
    DATA_MATRIX(Wpred);
    DATA_MATRIX(WxCLpred);
    //------------
    // PARAMETERS
    //------------
    PARAMETER_MATRIX(beta0);
    PARAMETER_VECTOR(betacond);
    PARAMETER_VECTOR(gammaw);
    PARAMETER_VECTOR(gammawcl);
    PARAMETER_VECTOR(a); // Cholesky factors
    PARAMETER_MATRIX(u0);
    //---------------------------
    // PRELIMINRARY CALCULATIONS
    //---------------------------
    // NUMBER OF OBSERVATIONS
    int n = Y.rows();
    // NUMBER OF CHOICES
    int m = Y.cols();
    // number of covariates
    int p = X.cols();
    // depth of the array
    int q = Xcond.cols(); // number of outermost dimensions - slabs
    // NUMBER OF PREDICTIONS
    int npred = Xpred.rows();
    // BETA WITH A COLUMN OF ZEROS FIRST
    matrix<Type> beta(p,m);
    for(int i = 0; i < p; i++){
      for(int j = 1; j < m; j++){ // leaves the first column as zeros
        beta(i,j) = beta0(i,j-1);
      }
    }
    // BULK WEIGHT EFFECTS
    matrix<Type> weight_effects(4,4);
    for(int i = 0; i < 4 ; i++){
      weight_effects(i,i) = gammaw(i);
    }
    // BULK WEIGHT x CL EFFECTS
    matrix<Type> weightcl_effects(4,4);
    for(int i = 0; i < 4 ; i++){
      weightcl_effects(i,i) = gammawcl(i);
    }
    // RANDOM EFFECTS WITH COLUMN OF ZEROS FIRST
    matrix<Type> u1(ngp,m);
    for(int i = 0; i < ngp; i++){
      for(int j = 1; j < m; j++){ // leaves the first column as zeros
        u1(i,j) = u0(i,j-1);
      }
    }
    // full matrix of random effects
    matrix<Type> u(n,m);
    for(int i = 0; i < n; i++){
        u.row(i) = u1.row(gp[i]);
    }
    // sigma
    matrix<Type> L(3,3);
    int ii = 0;
    for(int i = 0; i < 3; i++){
      for(int j = 0; j <= i; j++){
        L(i,j) = a(ii);
        ii += 1;
      }
    }
    matrix<Type> Sigma = L * L.transpose();
    //-----------
    // PROCEDURE
    //-----------
    matrix<Type> eta = X * beta + W * weight_effects + WxCL * weightcl_effects + u + Offset;
    // conditional variables
    for(int i = 0; i < q; i++){
      matrix<Type> etacondi = asMatrix(vector<Type>(Xcond.col(i) * betacond(i)), m, n); // transposed below 
        eta += etacondi.transpose();
    }
    // element-wise exponentiation
    matrix<Type> expeta(n,m);
    for(int i = 0; i < n; i++){
      for(int j = 0; j < m; j++){
          expeta(i,j) = exp(eta(i,j)); // can we do quicker?
      }
    }
    // row sums of exp(eta)
    matrix<Type> ones(m,1);
    ones.fill(1.0);
    matrix<Type> rowsumsexpeta =  expeta * ones;
    matrix<Type> P(n,m);
    for(int i = 0; i < n; i++){ 
      for(int j = 0; j < m; j++){ 
        P(i,j) = expeta(i,j) / rowsumsexpeta(i,0);
      }
    }
    Type nll = 0.0; // initialize negative log likelihood
    // random effects component
    for(int j = 0; j < ngp; j++){
      nll += MVNORM(Sigma, false)(vector<Type>(u0.row(j))); // note + as return negative log density
    }
    // observations component
    for(int i = 0; i < n; i++){ 
      nll -= dmultinom(vector<Type>(Y.row(i)), vector<Type>(P.row(i)), true);
    }
    // PREDICTION
    matrix<Type> etapred = Xpred * beta + Wpred * weight_effects + WxCLpred * weightcl_effects;
    ADREPORT(etapred);
    // conditional variables
    for(int i = 0; i < q; i++){
      matrix<Type> etacondpredi = asMatrix(vector<Type>(Xcondpred.col(i) * betacond(i)), m, npred); // transposed below 
      etapred += etacondpredi.transpose();
    }
   // ODDS
   matrix<Type> lodds_70_80 = etapred.col(0) - etapred.col(1);
   matrix<Type> p_80_70(npred,1);
   for(int i = 0; i < npred; i++){
     p_80_70(i,0) = exp(etapred(i,1)) / (exp(etapred(i,0)) + exp(etapred(i,1)));
   }
   matrix<Type> lodds_70_90 = etapred.col(0) - etapred.col(2);
   matrix<Type> lodds_70_100 = etapred.col(0) - etapred.col(3);
   matrix<Type> lodds_80_90 = etapred.col(1) - etapred.col(2);
   matrix<Type> lodds_80_100 = etapred.col(1) - etapred.col(3);
   matrix<Type> lodds_90_100 = etapred.col(2) - etapred.col(3);
   ADREPORT(lodds_70_80);
   ADREPORT(p_80_70);
   ADREPORT(lodds_70_90);
   ADREPORT(lodds_70_100);
   ADREPORT(lodds_80_90);
   ADREPORT(lodds_80_100);
   ADREPORT(lodds_90_100);
   return nll;
  }
"

write(tmb_multinomial_re, file = "multinomial_re_model.cpp")

compile("multinomial_re_model.cpp")

dyn.unload(dynlib("multinomial_re_model")) ## before loading to be safe
dyn.load(dynlib("multinomial_re_model"))

## DATA FOR TMB
## WEIGHTS
bulk.weights.mat <- as.matrix(our.lass.neph.cast[, c("m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")])
colnames(bulk.weights.mat) <- rownames(bulk.weights.mat) <- NULL
W <- scale.fun(bulk.weights.mat) ## to maintain relative weights of nets

## WEIGHTS x CL INTERACTION
WxCL <- scale.fun(bulk.weights.mat * our.lass.neph.cast[, "Carapace.length"])

## POSITION X CL INTERACTION
PIxCL <- PI * scale.fun(our.lass.neph.cast$Carapace.length)
SIxCL <- SI * scale.fun(our.lass.neph.cast$Carapace.length)
SOxCL <- SO * scale.fun(our.lass.neph.cast$Carapace.length)

## FILL THE CONDITIONAL ARRAY
Xcond.array <- array(NA, dim = c(n, m, 6))
Xcond.array[, , 1] <- PI
Xcond.array[, , 2] <- SI
Xcond.array[, , 3] <- SO
Xcond.array[, , 4] <- PIxCL
Xcond.array[, , 5] <- SIxCL
Xcond.array[, , 6] <- SOxCL

## PREDICTION OBJECTS FOR TMB
## X
npred <- 50
CLpred <- seq(min(our.lass.neph.cast$Carapace.length), max(our.lass.neph.cast$Carapace.length), length = npred)
CLpred.sc <- with(our.lass.neph.cast, (CLpred - mean(Carapace.length)) / sd(Carapace.length))
Xpred <- cbind(1, CLpred.sc)

## Xcond
Xcondpred <- array(NA, dim = c(npred, m, 6))

## proportions in given positions in the data
(our.lass.neph.cast[,net.names])

PIpred <- SIpred <- SOpred <- matrix(1/4, nr = npred, nc = 4) ## 1/4 as sums to one across the four levels
Xcondpred[, , 1] <- PIpred
Xcondpred[, , 2] <- PIpred
Xcondpred[, , 3] <- PIpred
Xcondpred[, , 4] <- PIpred * CLpred.sc
Xcondpred[, , 5] <- SIpred * CLpred.sc
Xcondpred[, , 6] <- SOpred * CLpred.sc

## W pred
W.names <- paste("m", seq(70,100, by = 10), "mm_Total.catch", sep = "")
W.means <- as.numeric(apply(our.lass.neph.cast[, W.names], 2, mean))
mean.W <- mean(unlist(our.lass.neph.cast[, W.names]))
sd.W <- sd(unlist(our.lass.neph.cast[, W.names]))
Wsc.means <- (W.means - mean.W) / sd.W
Wpred <- matrix(as.numeric(Wsc.means), nr = 1)[rep(1, npred), ]

## WxCL pred
mean.WxCL <- mean(unlist(our.lass.neph.cast[, W.names] * our.lass.neph.cast[, "Carapace.length"]))
sd.WxCL <- sd(unlist(our.lass.neph.cast[, W.names] * our.lass.neph.cast[, "Carapace.length"]))

WxCLpred <- (t(sapply(CLpred, FUN = function(x){ matrix(x * W.means)})) - mean.WxCL) / sd.WxCL

rm(obj)

##load("objects_for_Brian.RData")

obj <- MakeADFun(
    data = list(
        Y = Y,
        X = X,
        Xcond = Xcond.array,
        W = W,
        WxCL = WxCL,
        Offset = offset.mat,
        gp = gps - 1,
        ngp = ngp,
        ## prediction matrices
        Xpred = Xpred,
        Xcondpred = Xcondpred,
        Wpred = Wpred,
        WxCLpred = WxCLpred),
    parameters = list(
        beta0 = matrix(0, ncol = ncol(Y) - 1, nrow = ncol(X)),
        betacond = rep(0, dim(Xcond.array)[3]),
        gammaw = rep(0, 4),
        gammawcl = rep(0, 4),
        ##a = c(1,0,1,0,0,1), ## when not converging, sometimes need to change these values
        a = c(1.01, 0, 1.01, 0, 0, 1.01),
        u0 = matrix(0, ncol = ncol(Y) - 1, nrow = ngp)
        ),
    map = list(
        ## if need to fix any parameters
        beta0 = factor(c(1,2,3,NA,4,NA)),
        betacond  = factor(c(1,NA,2,3,NA,4)),
        gammaw = factor(c(1,2,3,4)),
        gammawcl = factor(c(1,2,NA,NA)),
        ##a = factor(rep(NA, 6)),
        ##u0 = factor(matrix(NA, ncol = ncol(Y) - 1, nrow = ngp))
        ),
    random = c("u0"),
    DLL = "multinomial_re_model"
    )

(opt<-do.call("optim", obj))
## (opt_nore<-do.call("optim", obj))

##(rep <- sdreport(obj))
(rep_nore<- sdreport(obj))

save(list = c("opt", "rep"), file = "best_fit.RData")
##save(list = c("opt_nore", "rep_nore"), file = "best_fit_nore.RData")


##----------------------
## Over-dispersion test
##----------------------
## likelihood obtained without random effects
ll <- -2105.376
## saturated likelihood - na.omit for log(0)
llsat.vec <- sapply(1:nrow(Y), function(z){
  lfactorial(sum(Y[z,])) - sum(lfactorial(Y[z,])) + sum(na.omit(Y[z,] * log(Y[z,] / sum(Y[z,]))))
  ##log(factorial(sum(Y[z,]))/prod(factorial(Y[z,])) * prod((Y[z,]/sum(Y[z,]))^Y[z,]))
})

llsat <- sum(llsat.vec)
## deviance
D <- -2 * (ll - llsat)

curve(dchisq(x, df = nrow(Y) * (ncol(Y) - 1) - length(opt$par)), from = 0, to = 3e3, col = "purple", lwd = 2)
abline(v = D, lty = 2)

pchisq(D, df = nrow(Y) * (ncol(Y) - 1) - length(opt$par), lower.tail = FALSE)



## simulation for distribution of the deviance
n <- 1e3
dev.vec <- rep(NA, n)

for(i in 1:n){
  print(i)
  X <- cbind(1, rnorm(400))
  beta <- cbind(0, matrix(rnorm(6), nrow = 2))
  eta <- X %*% beta
  p <- exp(eta) / rowSums(exp(eta))
  Y <- apply(p, 1, FUN = function(z){
    rmultinom(1, size = 100, prob = z)
  })
  Y <- t(Y)
  fit <- multinom(Y ~ X[,2], trace = FALSE)
  ##
  ll.vec <- sapply(1:nrow(Y), function(z){
    lfactorial(sum(Y[z,])) - sum(lfactorial(Y[z,])) + sum(Y[z,] * log(predict(fit, type = "prob")[z,]))
  })
  ll <- sum(ll.vec)
  ## saturated model - note na.omit
  llsat.vec <- sapply(1:nrow(Y), function(z){
    ##lfactorial(sum(Y[z,])) - sum(lfactorial(Y[z,])) + sum(na.omit(Y[z,] * log(Y[z,]/100)))
    log(factorial(sum(Y[z,]))/prod(factorial(Y[z,])) * prod((Y[z,]/100)^Y[z,]))
  })
  llsat <- sum(llsat.vec)
  ## library(MGLM)
  ## mglm.fit <- MGLMreg(Y ~ 1, dist="MN")
  dev.vec[i] <- -2 * (ll - llsat)
}

## note that changing the standard deviation of the parameters
## results in a closer fit to the predicted
hist(dev.vec, probability = TRUE, breaks = 30, border = "darkgrey", main = "")
curve(dchisq(x, df = nrow(Y) * (ncol(Y) - 1) - prod(dim(beta[,-1]))), add = TRUE, col = "purple", lwd = 2)

plot(density(dev.vec))


##-----------
## ADMB CODE
##-----------
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


## Conditional variable - catch per net
Xcond <- our.lass.neph.cast[, c("m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")]


## NOW SEND THESE OVER TO ADMB
n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)
q <- ncol(Xcond)
nchol <- (m-1)*(m)/2

beta0.start <- matrix(0, nrow = p, ncol = m - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)

## write the data out
datfile <- "../admbre/multinomialme.dat"
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
## testing
cat("# y \n", c(0.2,0.1) - c(1,-1), "\n", file = datfile, append = TRUE)
cat("# sigma \n", diag(2), "\n", file = datfile, append = TRUE)

## pinfile
datfile <- "../admbre/multinomialme.pin"
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

## run e vous ./multinomialme -l1 50000000 -l2 200000000 -l3 50000000

coef.admb <- read.table("../admbre/multinomialme.std", header = TRUE)

beta.admb <- matrix(coef.admb$value[coef.admb$name == "beta0"], ncol = ncol(X), byrow = FALSE)

betacond.admb <- coef.admb$value[coef.admb$name == "betacond"]

uhat.admb <- matrix(subset(coef.admb, name == "u0")$value, ncol = m-1, byrow = TRUE)

## Correlation of the random effects
a.admb <- subset(coef.admb, name == "a")$value

ii <- 1

L <- matrix(0, 3, 3)
for(i in 1:3){
    for(j in 1:i){
        L[i,j] <- a.admb[ii]
        ii <- ii + 1
    }
}

sigma.admb <- L %*% t(L)


## COMPARISON

## Compare
## beta
beta.admb; round(beta.tmb, 5)
plot(beta.admb, beta.tmb); abline(c(0,1))
## check filling correct
## plot(beta.admb, coef(mnom.fit)); abline(c(0,1))
## conditional beta
betacond.admb; betacond.tmb
## sigma 
sigma.admb; sigma.tmb
## random effects
plot(uhat.admb, uhat.tmb); abline(c(0,1))

## why do the likelihoods differ?

## use tmb final values as starting noeval values in admb
datfile <- "../admb/conlogit_re/multinomialme.pin"
cat("# beta0 \n", file = datfile)
write.table(t(beta.tmb), file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# betacond \n", betacond.tmb, "\n", file = datfile, append = TRUE)
##set.seed(10)
cat("# a \n", a.tmb, "\n", file = datfile, append = TRUE)
##cat("# a \n", rep(0.1, m-1), "\n", file = datfile, append = TRUE)
cat("# u0 \n", file = datfile, append = TRUE)
write.table(uhat.tmb, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## use tmb value in admb
obj2 <- MakeADFun(
    data = list(Y = Y, X = X, Xcond = Xcond.array, Offset = offset.mat, gp = gps - 1, ngp = ngp),
    parameters = list(
        beta0 = t(beta.admb),
        betacond = betacond.admb,
        a = a.admb,
        u0 = uhat.admb
        ),
    random = c("u0"),
    DLL = "multinomial_re_model")

obj2$fn()
## returns same likelihood as TMB so difference is in the way it is being calculated
## check the densities
dmultinom(c(1,2,3,4), p = rep(1,4)/4, log = TRUE)
library(mvtnorm)
dmvnorm(c(0.2,0.1), mean = c(1,-1), sigma = diag(2), log = TRUE)
