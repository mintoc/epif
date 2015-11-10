##-------------------
## Quad-rig plots
## CM, BB: Tue Sep 22 2015
## Note: bigger models run on the cluster
## CIs probably best estimated via:
## Sison, C.P and J. Glaz. 1995. Simultaneous confidence intervals and sample size determination for multinomial proportions. Journal of the American Statistical Association, 90:366-369. Paper available at http://tx.liberal.ntu.edu.tw/~purplewoo/Literature/!Methodology/!Distribution_SampleSize/SimultConfidIntervJASA.pdf
## No time to do that but reference
##-------------------
library(ggplot2)

## load the data
load("our_lass_data_objects.RData")

## load model objects
load("model_objects.RData")

## Get the proportions
count.mesh <- as.matrix(our.lass.neph.cast[, count.vars]/our.lass.neph.cast[, subsratio.vars])

prop.mesh <- prop.table(count.mesh, margin = 1)

m <- dim(prop.mesh)[1]

## make a dataframe of the proportions for ggplot
prop.mesh.df <- data.frame(
                  Mesh.Size = factor(rep(count.vars, each = m)),
                  Carapace.length = rep(our.lass.neph.cast$Carapace.length, times = 4),
                  fHAUL = rep(our.lass.neph.cast$fHAUL, times = 4),
                  Net.position = unlist(our.lass.neph.cast[, nc.vars]),
                  proportion = c(prop.mesh),
                  count = c(count.mesh))

prop.mesh.df$Mesh.Size <- factor(as.character(prop.mesh.df$Mesh.Size), levels = c("m70mm_Count", "m80mm_Count", "m90mm_Count", "m100mm_Count"))

prop.mesh.df$pch.cex <- log10(prop.mesh.df$count + 10)

levels(prop.mesh.df$Mesh.Size) <- c("70mm", "80mm", "90mm", "100mm")

prop.plot <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion, colour = fHAUL)) + 
  geom_point(aes(size = pch.cex)) + 
facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75)) + theme(legend.position = "none")

theme_set(theme_bw(base_size = 12))

png("../tex/figures/bubble_gum_plot_overall_v0.png", height = 7, width = 8, units = "in", res = 400)
prop.plot
dev.off()

## Plot to add overall fit back in
prop.plot.nocolour <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion)) + 
  ##geom_point(aes(size = pch.cex), colour = "#F8766D") +
  geom_point(aes(size = pch.cex), pch = 1, colour = "darkgrey") + 
facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75)) + theme(legend.position = "none")

## Plot to add by-haul fits in
## By Haul
prop.plot.haul <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion)) + 
  geom_point(aes(size = pch.cex), shape = 1, colour = "darkgrey") + 
  facet_grid(fHAUL ~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end")  + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75))
##scale_colour_gradientn(colours=rainbow(4))


##---------------
## INCLUDE FITS
##---------------
## read back in the results
coef.admb <- read.table("../admb/conlogit_re/conditional.std", header = TRUE)

m <- 4
uhat <- matrix(subset(coef.admb, name == "u0")$value, ncol = m-1, byrow = TRUE)

png("../tex/figures/multinomial_random_effects_v0.png", height = 6, width = 7, units = "in", res = 400)
par(mar = c(4, 4, 2, 2))
matplot(1:12, uhat, pch = c(2, 8, 19), col = 1, xlab = "", ylab = "", xaxt = "n")
axis(side = 1, at = 1:12)
abline(h = 0, lty = 2)
mtext(side = 1, line = 2.5, text = "Haul number")
mtext(side = 2, line = 2.5, text = "Log-odds random effects")
legend("topright", legend = c("80mm/70mm", "90mm/70mm", "100mm/70mm"), pch = c(2, 8, 19))
dev.off()

## Correlation of the random effects
a.hat <- subset(coef.admb, name == "a")$value

ii <- 1

L <- matrix(0, 3, 3)
for(i in 1:3){
  for(j in 1:i){
    L[i,j] <- a.hat[ii]
    ii <- ii + 1
  }
}

sigma.hat <- L %*% t(L)
cor(data.frame(uhat))

write.csv(x = round(cov2cor(sigma.hat), 2), "../data/correlation_u_hat.csv")

##
p <- dim(X)[2]
beta.hat <- cbind(0, matrix(subset(coef.admb, name == "beta0")$value, nrow = p, byrow = TRUE))
plot(coef(mnom.fit), t(beta.hat[,-1])); abline(c(0,1))

betacond.hat <- subset(coef.admb, name == "betacond")$value

cl.pred.vec <- Xpred[,"cl.vec"]

## predictions from admb
eta.hat.admb <- matrix(subset(coef.admb, name == "etapred")$value, ncol = 4, byrow = TRUE)
P.hat.admb <- exp(eta.hat.admb) / rowSums(exp(eta.hat.admb))
eta.se.admb <- matrix(subset(coef.admb, name == "etapred")$std.dev, ncol = 4, byrow = TRUE)

eta.lwr.admb <- eta.hat.admb + qnorm(0.025) * eta.se.admb
eta.upr.admb <- eta.hat.admb + qnorm(0.975) * eta.se.admb

P.admb.lwr <- exp(eta.lwr.admb) / rowSums(exp(eta.lwr.admb))
P.admb.upr <- exp(eta.upr.admb) / rowSums(exp(eta.upr.admb))

plot.pred.df.admb <- data.frame(
                       Mesh.Size = factor(rep(colnames(Y), each = nrow(Xpred))),
                       Carapace.length = rep(cl.pred.vec, times = 4),
                       proportion = c(P.hat.admb),
                       lwr = c(P.admb.lwr),
                       upr = c(P.admb.upr))

levels.switch <- c("m70mm_Count" = "70mm", "m80mm_Count" = "80mm", "m90mm_Count" = "90mm", "m100mm_Count" = "100mm")

levels(plot.pred.df.admb$Mesh.Size) <- levels.switch[levels(plot.pred.df.admb$Mesh.Size)]


png("../tex/figures/bubble_gum_plot_overall__with_fit_v0.png", height = 7, width = 8, units = "in", res = 400)
prop.plot.nocolour +
  ##geom_ribbon(data = plot.pred.df.admb, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.4) +
  geom_line(data = plot.pred.df.admb, aes(y = proportion), linetype = "solid") +
  geom_line(data = plot.pred.df.admb, aes(y = lwr), linetype = "dashed") +
  geom_line(data = plot.pred.df.admb, aes(y = upr), linetype = "dashed")
dev.off()

## predictions at the haul level
## note haul levels ordered even though obs unordered
## so reffs ordered
levels(our.lass.neph.cast$fHAUL)

Xcond.obs.df <- unique(our.lass.neph.cast[, c("fHAUL", "m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")])

Xcond.obs <- as.matrix(Xcond.obs.df[order(Xcond.obs.df$fHAUL), 2:5])

## first haul
eta.tmp <- apply(Xpred %*% beta.hat, 1, "+",
                 matrix(betacond.hat * Xcond.obs[1,] + c(0, uhat[1,]), nrow = 1)
                 ##matrix(betacond.hat * Xcond.obs[1,], nrow = 1)
                 )
eta.hat <- t(eta.tmp)
eta.df <- as.data.frame(eta.hat)
eta.df$fHAUL <- paste("H", 1, sep = "")

eta.df.all <- eta.df

for(i in 2:12){
  eta.tmp <- apply(Xpred %*% beta.hat, 1, "+",
                   matrix(betacond.hat * Xcond.obs[i,] + c(0, uhat[i,]), nrow = 1)
                   ##matrix(betacond.hat * Xcond.obs[i,], nrow = 1)
                   )
  eta.hat <- t(eta.tmp)
  eta.df <- as.data.frame(eta.hat)
  eta.df$fHAUL <- paste("H", i, sep = "")
  eta.df.all <- rbind(eta.df.all, eta.df)
}

P.all <- as.matrix(exp(eta.df.all[, 1:4]) / rowSums(exp(eta.df.all[, 1:4])))

fhaul.vec <- paste("H", 1:12, sep = "")

pred.df <- data.frame(proportion = c(P.all),
                      Carapace.length = rep(rep(cl.pred.vec, 12), 4),
                      fHAUL = rep(rep(fhaul.vec, each = length(cl.pred.vec)), 4),
                      Mesh.Size = rep(count.vars, each = length(cl.pred.vec) * 12))

levels(pred.df$Mesh.Size) <- levels.switch[levels(pred.df$Mesh.Size)]

##ggplot(pred.df, aes(x = Carapace.length, y = proportion)) + geom_line() + facet_grid(fHAUL ~ Mesh.Size)

png("../tex/figures/haul_level_predictions_v0.png", height = 10, width = 7, units = "in", res = 400)
prop.plot.haul + geom_line(data = pred.df) + theme(legend.position = "none")
dev.off()

## some of the haul fits look low?

