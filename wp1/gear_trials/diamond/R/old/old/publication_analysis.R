##-------------------
## Quad-rig analyses
## CM, BB: Tue Sep 22 2015
## Note: bigger models run on the cluster
##-------------------
library(gdata)

##----------
## OUR LASS 
##----------
our.lass.neph.dat <- read.xls("../data//2015 BIM Nephrops quad rig trials/Our Lass 2 70_80_90_100mm codends Irish Sea July 2015/Nephrops Raised Counts Our Lass 2 Irish Sea July 2015.xlsx", 
                              sheet = "All hauls",
                              stringsAsFactors = FALSE)

## removing haul 14 as only one tow with that net configuration
## check this again
our.lass.neph.dat <- subset(our.lass.neph.dat, Haul.No != 14)

## Each mesh size has a unique codend number
with(our.lass.neph.dat, table(Mesh.Size, Codend.No))

## position switched
with(our.lass.neph.dat, table(Mesh.Size, Net.position))

## Show the first 2 rows
head(our.lass.neph.dat, 2)

## Make the "HAUL" variable character
our.lass.neph.dat$HAUL <- paste("H", our.lass.neph.dat$Haul.No, sep ="")

## make some factor variables used in the analyses
our.lass.neph.dat$fHAUL <- factor(our.lass.neph.dat$HAUL, levels = unique(our.lass.neph.dat$HAUL))
our.lass.neph.dat$Mesh.Size <- factor(paste("m", our.lass.neph.dat$Mesh.Size, sep = ""))

## get count per length bin per haul by mesh size
## using the reshape package (makes it easier to process data)
library(reshape)

## variables to keep 
vars2keep <- c("Mesh.Size", "Carapace.length", "fHAUL", "Count")

## melt the data frame
our.lass.neph.melt <- melt(our.lass.neph.dat[, vars2keep], 
                  id = c("Mesh.Size", "Carapace.length", "fHAUL"))

## re-form the dataframe in required format 
our.lass.neph.cast <- cast(our.lass.neph.melt, Carapace.length + fHAUL ~ Mesh.Size  + variable)
our.lass.neph.cast <- our.lass.neph.cast[order(our.lass.neph.cast$fHAUL, our.lass.neph.cast$Carapace.length), ]
our.lass.neph.cast[is.na(our.lass.neph.cast)] <- 0

## check some
our.lass.neph.cast[our.lass.neph.cast$Carapace.length == 26 &
                   our.lass.neph.cast$fHAUL == "H10", "m70mm_Count"]

subset(our.lass.neph.dat,
       Mesh.Size == "m70mm" &
       Carapace.length == 26 &
       fHAUL == "H10")

## show the first few rows
head(our.lass.neph.cast, 2)

## format the subsampling ratio similarly

## unique raising factors per haul
rf.count <- with(our.lass.neph.dat, table(fHAUL, Overall.raising.factor, Mesh.Size))
apply(rf.count, 1, FUN = function(x){sum(x>0)})

## convert to sub-sampling ratio as in Celtic Warrior
our.lass.neph.dat$SUBSRATIO <- our.lass.neph.dat$Overall.Sampling.Ratio
vars2keep <- c("Mesh.Size", "fHAUL", "SUBSRATIO")
subs.melt <- melt(unique(our.lass.neph.dat[, vars2keep]), id = c("Mesh.Size", "fHAUL"))
subs.cast <- cast(subs.melt, fHAUL  ~ Mesh.Size + variable)

## check some
subs.cast[subs.cast$fHAUL == "H1", "m70mm_SUBSRATIO"]

subset(our.lass.neph.dat,
       Mesh.Size == "m70mm" &
       Carapace.length == 26 &
       fHAUL == "H1")

## get net position of each
vars2keep <- c("Mesh.Size", "fHAUL", "Net.position")
netpos.melt <- melt(unique(our.lass.neph.dat[, vars2keep]), id = c("Mesh.Size", "fHAUL"))
netpos.cast <- cast(netpos.melt, fHAUL  ~ Mesh.Size + variable)

## get bulk weight of each haul
vars2keep <- c("Mesh.Size", "fHAUL", "Total.catch")
bulk.melt <- melt(unique(our.lass.neph.dat[, vars2keep]), id = c("Mesh.Size", "fHAUL"))
bulk.cast <- cast(bulk.melt, fHAUL  ~ Mesh.Size + variable)


## merge counts and subsampling ratio back together 
our.lass.neph.cast0 <- merge(our.lass.neph.cast, subs.cast, by = "fHAUL", all.x = TRUE)

our.lass.neph.cast1 <- merge(our.lass.neph.cast0, bulk.cast, by = "fHAUL", all.x = TRUE)

our.lass.neph.cast <- merge(our.lass.neph.cast1, netpos.cast, by = "fHAUL", all.x = TRUE)

## show first few lines
head(our.lass.neph.cast, 2)

## Create "net configuration" variable
nc.vars <- c("m70mm_Net.position", "m80mm_Net.position", "m90mm_Net.position", "m100mm_Net.position")
our.lass.neph.cast$netconfig <- factor(paste("NC", apply(our.lass.neph.cast[, nc.vars], 1, paste, collapse = ""), sep =""))

## Extract the matrix of counts
count.vars <- c("m70mm_Count", "m80mm_Count", "m90mm_Count", "m100mm_Count")

neph.count.mat <- as.matrix(our.lass.neph.cast[, count.vars])

## Extract the matrix of subsampling ratios
subsratio.vars <- c("m70mm_SUBSRATIO", "m80mm_SUBSRATIO", "m90mm_SUBSRATIO", "m100mm_SUBSRATIO")

subsratio.mat <- as.matrix(our.lass.neph.cast[, subsratio.vars])

## Create the offset (NEED TO CHECK THIS)
offset.mat <- log(apply(subsratio.mat, 2, FUN = 
                        function(zz){zz/subsratio.mat[,1]}))

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

## MODELS

library(nnet)

## return to standardization with polynomials here 
## how about with quadratics?

polyfun <- poly(our.lass.neph.cast$Carapace.length, 2)

## using poly otherwise very strong correlations induced
mnom.fit <- multinom(neph.count.mat ~ 
                     Carapace.length +
                     netconfig +
                     poly(m70mm_Total.catch,2) +
                     poly(m80mm_Total.catch,2) +
                     poly(m90mm_Total.catch,2) +
                     poly(m100mm_Total.catch,2) +
                     offset(log(cbind(m70mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m80mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m90mm_SUBSRATIO/m70mm_SUBSRATIO,
                                      m100mm_SUBSRATIO/m70mm_SUBSRATIO))),
                     data = our.lass.neph.cast)

library(fields)
par.corr <- cov2cor(vcov(mnom.fit))
image.plot(par.corr)
hist(par.corr[lower.tri(par.corr)], breaks = 100, xlim = c(-1, 1), xaxs = "i")

poly70 <- poly(our.lass.neph.cast$m70mm_Total.catch, 2)
poly80 <- poly(our.lass.neph.cast$m80mm_Total.catch, 2)
poly90 <- poly(our.lass.neph.cast$m90mm_Total.catch, 2)
poly100 <- poly(our.lass.neph.cast$m100mm_Total.catch, 2)

par(mfrow = c(2,2), mar = c(2,2,1,1))
## 70mm
boxplot(resid(mnom.fit)[,1] ~ our.lass.neph.cast$fHAUL, notch = TRUE, col = "lightgrey"); abline(h = 0, lty =2)
legend("topleft", legend = "70mm")
## 80mm
boxplot(resid(mnom.fit)[,2] ~ our.lass.neph.cast$fHAUL, notch = TRUE, col = "lightgrey"); abline(h = 0, lty =2)
legend("topleft", legend = "80mm")
## 90mm
boxplot(resid(mnom.fit)[,3] ~ our.lass.neph.cast$fHAUL, notch = TRUE, col = "lightgrey"); abline(h = 0, lty =2)
legend("topleft", legend = "90mm")
## 100mm
boxplot(resid(mnom.fit)[,4] ~ our.lass.neph.cast$fHAUL, notch = TRUE, col = "lightgrey"); abline(h = 0, lty =2)
legend("topleft", legend = "100mm")

## ICC - not very strong for all covariates included
library(lme4)
fit <- lmer(resid(mnom.fit)[,4] ~ -1 + (1|our.lass.neph.cast$fHAUL), REML = FALSE)
sigmaa2 <- 0.1014^2
sigmae2 <- 0.1892^2
sigmaa2 / (sigmaa2 + sigmae2)

## ADMB-RE
Y <- neph.count.mat
X <- model.matrix(mnom.fit)
gps <- as.numeric(our.lass.neph.cast$fHAUL)
ngp <- length(unique(gps))

n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)
nchol <- (m-1)*(m)/2

beta0.start <- matrix(0, nrow = p, ncol = m - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)

## write the data out
datfile <- "../admb/multinomial_re/our_lass/multinomial.dat"
cat("# number of observations n \n", n, "\n", file = datfile)
cat("# number of categories m \n", m, "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", p, "\n", file = datfile, append = TRUE)
cat("# number of groups ngp \n", ngp, "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Y, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# groups \n", gps, "\n", file = datfile, append = TRUE)
cat("# Offset (log(q[i]/q[1])) \n", file = datfile, append = TRUE)
write.table(offset.mat, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## pinfile
datfile <- "../admb/multinomial_re/our_lass/multinomial.pin"
cat("# beta0 \n", file = datfile)
write.table(beta0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
##set.seed(10)
##cat("# a \n", rep(0.1, nchol), "\n", file = datfile, append = TRUE)
cat("# a \n", rep(0.1, m-1), "\n", file = datfile, append = TRUE)
cat("# u0 \n", file = datfile, append = TRUE)
write.table(u0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## read back in the results
coef.admb <- read.table("../admb/multinomial_re/our_lass/multinomial.std", header = TRUE)

uhat <- matrix(subset(coef.admb, name == "u0")$value, ncol = m-1, byrow = TRUE)

matplot(1:12, uhat)

plot(data.frame(uhat))

cor(data.frame(uhat))

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

beta.hat <- cbind(0, matrix(subset(coef.admb, name == "beta0")$value, nrow = p, byrow = TRUE))

plot(coef(mnom.fit), t(beta.hat[,-1])); abline(c(0,1))

## predictions
cl.pred.vec <- seq(min(our.lass.neph.cast$Carapace.length), max(our.lass.neph.cast$Carapace.length), length = 100)

q <- length(cl.pred.vec)

mean.bulk <- mean(unlist(bulk.cast[, -1]))
##mean.bulk <- 500
mean.bulk.vec <- apply(bulk.cast[,-1][,c("m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")], 2, mean)
##mean.bulk <- 0
## check the 1/4
levels(as.factor(our.lass.neph.cast$netconfig))

X.pred <- cbind(1, cl.pred.vec, 1/3, 1/3, 1/3,
                predict(poly70, mean.bulk)[1], predict(poly70, mean.bulk)[2],
                predict(poly80, mean.bulk)[1], predict(poly80, mean.bulk)[2],
                predict(poly90, mean.bulk)[1], predict(poly90, mean.bulk)[2],
                predict(poly100, mean.bulk)[1], predict(poly100, mean.bulk)[2]
                )
##X.pred <- cbind(1, cl.pred.vec, mean.bulk.vec[1], mean.bulk.vec[2], mean.bulk.vec[3], mean.bulk.vec[4])
##X.pred <- cbind(1, cl.pred.vec)
##X.pred <- cbind(1, 1/3, 1/3, 1/3, mean.bulk, mean.bulk, mean.bulk, mean.bulk)
## X.pred <- X.pred[rep(1,length(cl.pred.vec)),]

eta.hat <- X.pred %*% beta.hat

P.hat <- exp(eta.hat) / rowSums(exp(eta.hat))

## predictions from multinom change with each prediction - not good

eta.hat.mnom <- X.pred %*% cbind(0, t(coef(mnom.fit)))

P.hat.mnom <- exp(eta.hat.mnom) / rowSums(exp(eta.hat.mnom))

plot.pred.df <- data.frame(
                  Mesh.Size = factor(rep(colnames(Y), each = q)),
                  Carapace.length = rep(cl.pred.vec, times = 4),
                  proportion = c(P.hat),
                  proportion.mnom = c(P.hat.mnom))

prop.plot.nohaul + geom_line(data = plot.pred.df) + geom_line(data = plot.pred.df, aes(y = proportion.mnom), colour = "orange", linetype = "dashed")

## with colours for haul
eta.hat <- t(apply(X.pred %*% beta.hat, 1, "+", c(0, uhat[1,])))

for(i in 2:12){
  eta.hat <- rbind(eta.hat, t(apply(X.pred %*% beta.hat, 1, "+", c(0, uhat[i,]))))
}

P.hat <- exp(eta.hat) / rowSums(exp(eta.hat))

plot.pred.df <- data.frame(
                  Mesh.Size = factor(rep(colnames(Y), each = q * 12)),
                  Carapace.length = rep(cl.pred.vec, times = 4 * 12),
                  proportion = c(P.hat),
                  eta = c(eta.hat),
                  fHAUL = rep(rep(paste("H", 1:12, sep = ""), each = q),m))
##lower = c(pred.lower),
##upper = c(pred.upper))

## show parallel on linear predictor scale
ggplot(plot.pred.df, aes(x = Carapace.length, y = eta, group = fHAUL)) + geom_line() + facet_wrap(~ Mesh.Size)

pdf("../tex/figures/bubble_gum_byhaul_with_fits_plots.pdf", height = 10, width = 9)
prop.plot.haul + geom_line(data = plot.pred.df) ##+ geom_line(data = plot.pred.df, aes(y = proportion.mnom), colour = "orange", linetype = "dashed")
dev.off()

prop.plot + geom_line(data = plot.pred.df) ##+ geom_line(data = plot.pred.df, aes(y = proportion.mnom), colour = "orange", linetype = "dashed")

ggplot(plot.pred.df, aes(x = Carapace.length, y = proportion, colour = fHAUL)) + 
  geom_line() + facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") 

## generate some random effects to see if crossover can occur

x <- seq(0,10)

u <- cbind(0, mvrnorm(10, c(0,0,0), sigma.hat))
(eta.hat <- t(apply(cbind(1, x) %*% beta.hat, 1, "+", u[1,])))
P.hat <- exp(eta.hat) / (rowSums(exp(eta.hat)))



##---------
## SANDBOX 
##---------
## TEST THE SUB-SAMPLING RATIO
## MAKE UP SOME SUB-SAMPLED DATA
Y.samp <- t(rmultinom(10, size = 40, prob = c(0.1, 0.2, 0.4, 0.3)))

p.samp <- colSums(Y.samp) / sum(Y.samp)

n <- nrow(Y.samp)
m <- ncol(Y)
subsratio <- matrix(NA, nrow = n, ncol = m)
subsratio[,1] <- 0.1
subsratio[,2] <- 0.2
subsratio[,3] <- 0.5
subsratio[,4] <- 0.1
offset.mat <- log(apply(subsratio, 2, FUN = 
                        function(zz){zz/subsratio[,1]}))

Y.true <- Y.samp / subsratio
p.true <- colSums(Y.true) / sum(Y.true)
mnom.fit <- multinom(Y.samp ~ 1 + offset(offset.mat))
test <- exp(matrix(c(0, coef(mnom.fit)[,1]), ncol = 1))/sum(exp(matrix(c(0, coef(mnom.fit)[,1]), ncol = 1)))

## ADMB-RE
##Y <- neph.count.mat
Y <- Y.samp
X <- model.matrix(mnom.fit)
gps <- as.numeric(our.lass.neph.cast$fHAUL)

n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)
nchol <- (m-1)*(m)/2

beta0.start <- matrix(0, nrow = p, ncol = m - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)

gps <- rep(1, n)
ngp <- length(unique(gps))

## write the data out
datfile <- "../admb/multinomial_re/our_lass/multinomial.dat"
cat("# number of observations n \n", n, "\n", file = datfile)
cat("# number of categories m \n", m, "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", p, "\n", file = datfile, append = TRUE)
cat("# number of groups ngp \n", ngp, "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Y, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
cat("# groups \n", gps, "\n", file = datfile, append = TRUE)
cat("# Offset (log(q[i]/q[1])) \n", file = datfile, append = TRUE)
write.table(offset.mat, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## pinfile
datfile <- "../admb/multinomial_re/our_lass/multinomial.pin"
cat("# beta0 \n", file = datfile)
write.table(beta0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)
##set.seed(10)
##cat("# a \n", rnorm(nchol), "\n", file = datfile, append = TRUE)
cat("# a \n", rep(0.1, nchol), "\n", file = datfile, append = TRUE)
cat("# u0 \n", file = datfile, append = TRUE)
write.table(u0.start, file = datfile, append = TRUE, col.names = FALSE, row.names = FALSE)

## read back in the results
coef.admb <- read.table("../admb/multinomial_re/our_lass/multinomial.std", header = TRUE)

uhat <- matrix(subset(coef.admb, name == "u0")$value, ncol = ngp, byrow = TRUE)

matplot(1:13, uhat)

beta.hat <- cbind(0, matrix(subset(coef.admb, name == "beta0")$value, nrow = p, byrow = TRUE))

exp(beta.hat) / sum(exp(beta.hat))


