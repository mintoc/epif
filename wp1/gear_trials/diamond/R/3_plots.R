##-------------------
## Quad-rig plots
## CM, BB: Tue Sep 22 2015
## Note: bigger models run on the cluster
##-------------------



## PLOT

library(ggplot2)

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

prop.plot <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion, colour = fHAUL)) + 
  geom_point(aes(size = pch.cex)) + 
facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75))
##scale_colour_gradientn(colours=rainbow(4))

prop.plot.nohaul <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion)) + 
  geom_point(aes(size = pch.cex), colour = "#F8766D") + 
facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75))


theme_set(theme_bw(base_size = 16))

pdf("../tex/figures/bubble_gum_plots.pdf", height = 7, width = 9)
print(prop.plot)
dev.off()

prop.plot.haul <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion, colour = fHAUL)) + 
  geom_point(aes(size = pch.cex)) + 
  facet_grid(fHAUL ~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end")  + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75))
##scale_colour_gradientn(colours=rainbow(4))

pdf("../tex/figures/bubble_gum_byhaul_plots.pdf", height = 10, width = 9)
prop.plot.haul
dev.off()

##---------------
## FITS ADDED IN
##---------------
## read back in the results
coef.admb <- read.table("../admb/conlogit_re/conditional.std", header = TRUE)

uhat <- matrix(subset(coef.admb, name == "u0")$value, ncol = m-1, byrow = TRUE)

matplot(1:12, uhat)

plot(data.frame(uhat))

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

cov2cor(sigma.hat)
cor(data.frame(uhat))

beta.hat <- cbind(0, matrix(subset(coef.admb, name == "beta0")$value, nrow = p, byrow = TRUE))

betacond.hat <- subset(coef.admb, name == "betacond")$value

plot(coef(mnom.fit), t(beta.hat[,-1])); abline(c(0,1))


## fixed effect predictions
cl.pred.vec <- seq(min(our.lass.neph.cast$Carapace.length), max(our.lass.neph.cast$Carapace.length), length = 10)

##q <- length(cl.pred.vec)

mean.bulk <- 373

X.pred <- cbind(1, cl.pred.vec, 1/3, 1/3, 1/3)

Xcond.pred <- matrix(mean.bulk, nrow = nrow(X.pred), ncol = q)

eta.hat <- X.pred %*% beta.hat + betacond.hat * Xcond.pred

P.hat <- exp(eta.hat) / rowSums(exp(eta.hat))

## predictions from multinom change with each prediction - not good

eta.hat.mnom <- X.pred %*% cbind(0, t(coef(mnom.fit)))

P.hat.mnom <- exp(eta.hat.mnom) / rowSums(exp(eta.hat.mnom))

plot.pred.df <- data.frame(
                  Mesh.Size = factor(rep(colnames(Y), each = nrow(X.pred))),
                  Carapace.length = rep(cl.pred.vec, times = 4),
                  proportion = c(P.hat),
                  proportion.mnom = c(P.hat.mnom))

prop.plot.nohaul +
  geom_line(data = plot.pred.df, aes(x = Carapace.length, y = proportion)) ##+
  ##geom_line(data = plot.pred.df, aes(y = proportion.mnom), colour = "black", linetype = "dashed")

## predictions at the haul level
## note haul levels ordered even though obs unordered
## so reffs ordered
levels(our.lass.neph.cast$fHAUL)

Xcond.obs.df <- unique(our.lass.neph.cast[, c("fHAUL", "m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")])

Xcond.obs <- Xcond.obs.df[order(Xcond.obs.df$fHAUL), 2:5]



eta.hat <- X.pred %*% beta.hat + betacond.hat * Xcond.obs


betacond.hat * Xcond.obs + cbind(0, uhat)

eta.hat <- X.pred %*% beta.hat + betacond.hat * Xcond.obs



