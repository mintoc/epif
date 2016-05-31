

## poisson-multinomial link
##miners.glm5 <- glm(count~period+resp,data=miners.long,family=poisson)
miners.glm5 <- glm(count ~ -1 + resp, data = miners.long, family = poisson)

## same model in multinom
library(nnet)

Y <- with(miners, cbind(n, m, s))

mnom0 <- multinom(Y ~ 1, data = miners)

phat <- predict(mnom0, type = "prob")[1,]

mnom.ll <- sum(apply(Y, 1, FUN = function(x){dmultinom(x, size = sum(x), prob = phat, log = TRUE)}))

pred.rate <- predict(miners.glm5, newdata = data.frame(resp = c("n", "m", "s")))

prop.table(apply(miners[, c("n", "m", "s")], 2, sum))

mu.hat <- exp(coef(miners.glm5))

mu.hat / sum(mu.hat)

phi.hat <- sum(mu.hat)

n <- sum(Y)

- phi.hat + n * log(phi.hat) - lfactorial(n)

mnom.ll + logLik(miners.glm5)


mui.hat <- predict(miners.glm5, type = "response")

sum(- mui.hat + miners.long$count * log(mui.hat) - lfactorial(miners.long$count))


n <- rowSums(Y)

Prn <- sum(-phi.hat + n * log(phi.hat) - lfactorial(n))

logLik(miners.glm5) - Prn


