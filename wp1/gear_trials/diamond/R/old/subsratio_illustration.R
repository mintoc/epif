##------------------------------------
## Illustration of sub-sampling ratio
## CM:
##
##------------------------------------

## true counts
n <- c(10, 40, 30, 20)

## true proportions
p <- n / sum(n)

## subsampling ratios
subs <- c(.1, 1, 1, 1)

## observed counts 
m <- n * subs

m / sum(m)

(m / subs) / sum((m / subs))

## what it should be 
library(nnet)
fit0 <- multinom(t(matrix(n)) ~ 1, reltol = 1.0e-16)
predict(fit0, type = "probs")
## same as
exp(rbind(0, coef(fit0))) / sum(exp(rbind(0, coef(fit0))))

## what it looks to be
fit1 <- multinom(t(matrix(m)) ~ 1, reltol = 1.0e-16)
predict(fit1, type = "probs")
## same as
exp(rbind(0, coef(fit1))) / sum(exp(rbind(0, coef(fit1))))

## with an offset
offset.mat <- matrix(log(subs / subs[1]), ncol = 4)
## what it looks to be
fit2 <- multinom(t(matrix(m)) ~ 1 + offset(offset.mat), reltol = 1.0e-16, )
predict(fit2, type = "probs") ## not right
## coefs right though
exp(rbind(0, coef(fit2))) / sum(exp(rbind(0, coef(fit2))))



