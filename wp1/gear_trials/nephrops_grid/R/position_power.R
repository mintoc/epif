##------------------------------------------------------
## Preliminary power analysis of position in grid trial
## CM: Mon Jul 13 2015
## 
##------------------------------------------------------
## background
## Simulation to see which configuration might give the best power to detect differences for the grid trial. Design settings to be tested:
##   - Effects of the number of times the position of the nets is switched.
##   - Order of the switching, e.g., if you were to make one switch, what would it be to give the greatest contrast.

## There are a large number of factors that could influence these (combinations could be huge) but the idea is to keep it relatively simple, based on the data from other trials (to get an idea of the various variances) and then provide feedback on a limited number of scenarios ye envisage (change of grounds, logistics of changing the nets, etc.).

##----------
## SIMULATE
##----------
## True population - potentially ground-specific
## True selectivity of the gears
## Between-haul variability
## Within-haul variability
## Individual net effects
## Net position effects


## TRUE POPULATION

## median carapace lengths
mu.l <- c(10, 15, 20, 25)

## cv of carapace length at age 
cv.l <- 0.3

## proportion of the population composed of given component/age
p <- c(0.4, 0.3, 0.2, 0.1)
 
## population density at age
length.vec <- seq(0, 60, length = 1e3)

d1 <- p[1] * dlnorm(length.vec, log(mu.l[1]), sdlog = cv.l)
d2 <- p[2] * dlnorm(length.vec, log(mu.l[2]), sdlog = cv.l)
d3 <- p[3] * dlnorm(length.vec, log(mu.l[3]), sdlog = cv.l)
d4 <- p[4] * dlnorm(length.vec, log(mu.l[4]), sdlog = cv.l)

d <- rbind(d1, d2, d3, d4)

plot(length.vec, apply(d, 2, sum), type = "n", xlab = "", ylab = "")
matlines(length.vec, t(d), col = "purple", lty = 1)
lines(length.vec, apply(d, 2, sum), col = "forestgreen", lwd = 1.5)
legend("topright", legend = c("Component (age) distribution", "Overall distribution"), lty = 1, lwd = c(1, 1.5), col = c("purple", "forestgreen"), bty = "n")

curve(, from = 0, to = 60, n = 1e3)
curve(p[2] * dlnorm(length.vec, log(mu.l[2]), sdlog = cv.l), add = TRUE, n = 1e3)
curve(p[3] * dlnorm(length.vec, log(mu.l[3]), sdlog = cv.l), add = TRUE, n = 1e3)
curve(p[4] * dlnorm(length.vec, log(mu.l[4]), sdlog = cv.l), add = TRUE, n = 1e3)


