

a <- 100
b <- -0.1



p <- 

q <- 1:100
p0 <- 5
dp.dq <- -0.05

p <- p0 + dp.dq * q

plot(q, p, type = "l")

dq.dp <- 1 / dp.dq

e <- abs(dq.dp * (p / q))

par(mar = c(4, 4, 2, 4))
plot(q, p, type = "l", xlab = "Quantity", ylab = "Price")
par(new = TRUE)
plot(q, e, type = "l", lty = 2, axes = FALSE, ylab = "", xlab = "")
axis(side = 4)
mtext(side = 4, line = 2.5, text = "Elasticity")
legend("topright", legend = c("Price", "Elasticity"), lty = 1:2, bty = "n")
abline(h = 1, col = "grey")

