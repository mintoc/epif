##----------------------
## Attempt Poisson approx to multinomial with random effects
## CM: Fri Jul 31 2015
## Notes: compare random effects with just two categories - is it same as glmer for binomial
##----------------------

## simulate some repeated measures multinomial data

ntrials <- 50
nobs.per.trial <- 10

## covariate matrix
X <- cbind(1, rnorm(ntrials))

## parameters for three non-baseline categories
Beta <- matrix(c(3, 1, 2, 0.2, 0, -0.3), nrow = 3, ncol = 2)
Denominator <-  1 + exp(X %*% Beta[1,]) + exp(X %*% Beta[2,]) + exp(X %*% Beta[3,])

## Probabilities
Prob <- cbind(1/Denominator, exp(X %*% Beta[1,])/Denominator, exp(X %*% Beta[2,])/Denominator, exp(X %*% Beta[3,])/Denominator)

X.order <- X[order(X[,2]),2]
Prob.order <- Prob[order(X[,2]), ]

matplot(X.order, Prob.order, type = "l", lty = 1)

## Simulate some data
Y <- t(apply(Prob, 1, rmultinom, n = 1, size = 100))

matpoints(X[, 2], prop.table(Y, 1), pch = 1)

## Fit in nnet
library(nnet)

X.obs <- X[, 2]

m1 <- multinom(Y ~ X.obs)

coef(m1) ; Beta

## overlay predictions
X.pred <- seq(min(X[, 2]), max(X[, 2]), length = 100)

p1 <- predict(m1, newdata = data.frame(X.obs = X.pred), type = "probs")

matlines(X.pred, p1, type = "l", lty = 2)

## Fit using Poisson
## convert the data
long.dat <- data.frame(Y = c(Y),
                       category = rep(paste("cat", 1:4, sep = ""), each = ntrials),
                       ##obs = factor(rep(1:ntrials, each = 4)),
                       X = rep(X[, 2], 4))

pois1 <- glm(Y ~ category * X, family = poisson, data = long.dat)

coef(m1) ; coef(pois1); Beta

## predictions from this model
N.pred <- cbind(predict(pois1, newdata = data.frame(X = X.pred, category = "cat1"), type = "response"),
                predict(pois1, newdata = data.frame(X = X.pred, category = "cat2"), type = "response"),
                predict(pois1, newdata = data.frame(X = X.pred, category = "cat3"), type = "response"),
                predict(pois1, newdata = data.frame(X = X.pred, category = "cat4"), type = "response"))

p2 <- prop.table(N.pred, 1)
matlines(X.pred, p2, type = "l", lty = 4)

pois1 <- glmer(Y ~ category * X, random = ~ category | Haul,  family = poisson, data = long.dat)


##---------
## SANDBOX 
##---------


## car preference data
carpref = read.csv("CarPref.csv")
print(carpref)
print(summary(carpref))

# Create a multinomial response observation factor
# We can see that consecutive batches of four counts consist of one
# multinomial observation.

carpref$ResponseObs = factor(rep(1:6,rep(4,6)))
print(carpref)

# by default, "Important" is the reference level in the treatment contrast 
# coding because it is first in alphabetical order, and this is a bad choice.
# Force Not.Important to be the reference level

carpref$Color = relevel(carpref$Color, ref="Not.Important")
print(summary(carpref))  # much better!

# fit the multinomial model using the Poisson "trick"

# model 1 has Degree and Race appear additively
carpreffit1 = glm(Y ~ -1 + ResponseObs + Color*(Gender + Size),
   family=poisson, data=carpref)
print(summary(carpreffit1))

# model 2 includes only Race
carpreffit2 = glm(Y ~ -1 + ResponseObs + Color*(Size),
   family=poisson, data=carpref)
print(summary(carpreffit2))

# compare model fits
print(anova(carpreffit2,carpreffit1, test="Chi"))

# use function "multinom" in nnet library to fit 1st model
# Data needs to be organized differently though - see carprefexp.dat below

library(nnet)

carpref.exp = read.csv("CarPrefExp.csv")
print(carpref.exp) # each multinomial observation on its own line
carpref.exp.fit1 = multinom(
  cbind(Not.Important,Little.Importance,Important,Very.Important) ~
     Gender + Size, carpref.exp)
summary(carpref.exp.fit1)

# estimates are based on a neural network algorithm.
# they're not exactly MLEs!

# example of prediction: For "Male" and "Large"
predict(carpref.exp.fit1,
  newdata=data.frame(Gender="Male",Size="Large"),type="probs")

rm(carpref,carpreffit1,carpreffit2,carpref.exp,carpref.exp.fit1)

