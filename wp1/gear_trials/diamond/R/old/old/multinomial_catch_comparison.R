##------------------------------------
## Try a multinomial catch comparison
## CM: Wed May 27 2015
##
##------------------------------------

library(gdata)

setwd("../data")

neph.dat <- read.xls("Celtic Warrior Diamond mesh July 2014 Celtic Sea.xls", 
                     sheet = "Nephrops Lengths",
                     stringsAsFactors = FALSE)

## remove Haul 22, as no recordings for 90mm 
neph.dat <- subset(neph.dat, HAUL != 22)

## Show the first 2 rows
head(neph.dat, 2)

## Change the carapace length name
names(neph.dat)[names(neph.dat) == "Carapace.Length..mm.."] <- "Carapace.Length"

## Make the "HAUL" variable character
neph.dat$HAUL <- paste("H", neph.dat$HAUL, sep ="")

## make some factor variables variable
neph.dat$fHAUL <- factor(neph.dat$HAUL, levels = unique(neph.dat$HAUL))
neph.dat$fMesh.Size <- factor(neph.dat$Mesh.Size, levels = unique(neph.dat$Mesh.Size))

##-------------
## MULTINOMIAL 
##-------------
library(ggplot2)
library(nnet)
library(reshape)

## get count per length bin per haul per mesh size
vars2keep <- c("fMesh.Size", "Carapace.Length", "fHAUL", "COUNT")

neph.melt <- melt(neph.dat[, vars2keep], id = c("fMesh.Size", "Carapace.Length", "fHAUL"))
neph.cast <- cast(neph.melt, Carapace.Length + fHAUL ~ fMesh.Size  + variable)
neph.cast <- neph.cast[order(neph.cast$fHAUL, neph.cast$Carapace.Length), ]
neph.cast[is.na(neph.cast)] <- 0

## get the SUBSRATIO
vars2keep <- c("fMesh.Size", "fHAUL", "SUBSRATIO")
subs.melt <- melt(unique(neph.dat[, vars2keep]), id = c("fMesh.Size", "fHAUL"))
subs.cast <- cast(subs.melt, fHAUL  ~ fMesh.Size + variable)

## merge back in with nephrops counts
neph.cast <- merge(neph.cast, subs.cast, by = "fHAUL", all.x = TRUE)

## get proportion of total per mesh size
count.mesh <- as.matrix(neph.cast[, c("70mm_COUNT", "80mm_COUNT", "90mm_COUNT", "100mm_COUNT")])
prop.mesh <- prop.table(count.mesh, margin = 1)

m <- dim(prop.mesh)[1]

prop.mesh.df <- data.frame(Mesh.Size = factor(rep(c("70mm", "80mm", "90mm", "100mm"), each = m), levels = c("70mm", "80mm", "90mm", "100mm")),
                           Carapace.Length = rep(neph.cast$Carapace.Length, times = 4),
                           proportion = c(prop.mesh),
                           count = c(count.mesh))

library(ggplot2)

ggplot(prop.mesh.df, aes(x = Carapace.Length, y = proportion)) + geom_point(colour = "darkblue", alpha = 0.1, aes(size = log(count))) + facet_wrap(~ Mesh.Size) + ylab("Proportion of Nephrops per cod-end")

## simple multinomial for the proportions
count.vars <- c("70mm_COUNT", "80mm_COUNT", "90mm_COUNT", "100mm_COUNT")
neph.count.mat <- as.matrix(neph.cast[, count.vars])
colnames(neph.count.mat) <- c("70mm_COUNT", "80mm_COUNT", "90mm_COUNT", "100mm_COUNT")
## first fit
subsratio.mat <- as.matrix(neph.cast[, c("70mm_SUBSRATIO", "80mm_SUBSRATIO", "90mm_SUBSRATIO", "100mm_SUBSRATIO")])

offset.mat <- log(apply(subsratio.mat, 2, FUN = function(zz){zz/subsratio.mat[,1]}))

mnom0 <- multinom(neph.count.mat ~ 1 + offset(offset.mat))

## include carapace length third order polynomial (based on AIC and BIC)
## first scale it to range between zero and one
max.length <- max(neph.cast$Carapace.Length)
neph.cast$Carapace.Length <- neph.cast$Carapace.Length/max.length
neph.cast$Carapace.Length2 <- neph.cast$Carapace.Length^2
neph.cast$Carapace.Length3 <- neph.cast$Carapace.Length^3
## 
mnom.length <- multinom(neph.count.mat ~ Carapace.Length + Carapace.Length2 + Carapace.Length3 + offset(offset.mat), data = neph.cast)
##mnom.length <- multinom(neph.count.mat ~ Carapace.Length + Carapace.Length2 + Carapace.Length3, data = neph.cast)

## get predictions manually
## CIs not defined in multinomial context but let's try
beta.mu <- c(t(coef(mnom.length)))
Sigma <- vcov(mnom.length)

pred.length <- seq(min(neph.cast$Carapace.Length), max(neph.cast$Carapace.Length), length = 100)

nlength <- 100
nresamp <- 100
pred.array <- array(NA, dim = c(nlength, 4, nresamp))

library(mvtnorm)

X <- cbind(1, pred.length, pred.length^2, pred.length^3)

for(i in 1:nresamp){
  print(i)
  beta <- matrix(rmvnorm(1, mean = beta.mu, sigma = Sigma), nrow = 3, byrow = TRUE)
  p80 <- exp(X %*% matrix(beta[1,]))/(1 + rowSums(exp(X %*% t(beta))))
  p90 <- exp(X %*% matrix(beta[2,]))/(1 + rowSums(exp(X %*% t(beta))))
  p100 <- exp(X %*% matrix(beta[3,]))/(1 + rowSums(exp(X %*% t(beta))))
  p70 <- 1 - p80 - p90 - p100
  pred.p <- cbind(p70, p80, p90, p100)
  pred.array[ , , i] <- pred.p
  rm(pred.p)
}

pred.mu <- apply(pred.array, c(1, 2), mean)
pred.upper <- apply(pred.array, c(1, 2), quantile, p = 0.975)
pred.lower <- apply(pred.array, c(1, 2), quantile, p = 0.025)

m <- dim(pred.mu)[1]

pred.ci.df <- data.frame(
               Mesh.Size = factor(rep(c("70mm", "80mm", "90mm", "100mm"), each = m), levels = c("70mm", "80mm", "90mm", "100mm")),
               Carapace.Length = rep(pred.length * max.length, times = 4),
               proportion = c(pred.mu),
               lower = c(pred.lower),
               upper = c(pred.upper))

p <- ggplot(prop.mesh.df, aes(x = Carapace.Length, y = proportion)) + geom_point(colour = "#F8766D", alpha = 0.2, aes(size = log(count))) + facet_wrap(~ Mesh.Size) + ylab("Proportion of Nephrops per cod-end")

p + geom_ribbon(data=pred.ci.df, aes(ymin = lower, ymax = upper), alpha=0.3, fill = "blue") + geom_line(data = pred.ci.df, aes(x = Carapace.Length, y = proportion), col = "navy", size = 0.5) + geom_hline(aes(yintercept = 0.25), linetype = "dashed")
