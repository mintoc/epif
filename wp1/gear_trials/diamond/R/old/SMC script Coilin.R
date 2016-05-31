##--------------------------
## Multinomial code without random effects
## CM, DB: Feb 8th 2016
##--------------------------

library(gdata)

## Nephrops length data
##smc.dat <- read.csv("C:/Users/browned/Desktop/R/Data/StellaNovaSMC.csv", header=TRUE, stringsAsFactors=FALSE)
smc.dat <- read.csv("StellaNovaSMC.csv", header=TRUE, stringsAsFactors=FALSE)

## Make the "HAUL" variable character
smc.dat$HAUL <- paste("H", smc.dat$HAUL, sep ="")

## make some factor variables used in the analyses
smc.dat$fHAUL <- factor(smc.dat$HAUL, levels = unique(smc.dat$HAUL))
smc.dat$COMPARTMENT <- factor(smc.dat$COMPARTMENT)

## change a few names - this will work even if the name changes from the original
## if in doubt about names use names(grid.neph.dat) to find out what they are
names(smc.dat)[names(smc.dat) == "CARAPACE.LENGTH"] <- "Carapace.length"
names(smc.dat)[names(smc.dat) == "BULK.CATCH"] <- "Bulk.weight"

## Show the first 2 rows
head(smc.dat, 2)

## Bring in rotation data (used further down)
rotation.dat <-
  read.csv(
    "C:/Users/browned/Desktop/R/Data/Rotations.csv", 
    stringsAsFactors = FALSE, nrows = 13)

rotation.dat <-
  read.csv(
    "Rotations.csv", 
    stringsAsFactors = FALSE, nrows = 13)


##rotation.dat <- subset(rotation.dat, Haul.. != 4)

## need to create a rotation variable
## use short codes for rotation names
## starboard outside - added extra column in virtual spreadsheet
rotation.dat$SO <- NA
## now let's fill in the short names
rotation.dat$SO[rotation.dat$Starboard.Outside == "SMC45"] <- "SMC45"
rotation.dat$SO[rotation.dat$Starboard.Outside == "SMC55"] <- "SMC55"
rotation.dat$SO[rotation.dat$Starboard.Outside == "SMC65"] <- "SMC55"
rotation.dat$SO[rotation.dat$Starboard.Outside == "DMC75"] <- "DMC75"

## starboard inside
rotation.dat$SI <- NA
rotation.dat$SI[rotation.dat$Starboard.Inside == "SMC45"] <- "SMC45"
rotation.dat$SI[rotation.dat$Starboard.Inside == "SMC55"] <- "SMC55"
rotation.dat$SI[rotation.dat$Starboard.Inside == "SMC65"] <- "SMC55"
rotation.dat$SI[rotation.dat$Starboard.Inside == "DMC75"] <- "DMC75"

## Port outside
rotation.dat$PO <- NA
rotation.dat$PO[rotation.dat$Port.Outside == "SMC45"] <- "SMC45"
rotation.dat$PO[rotation.dat$Port.Outside == "SMC55"] <- "SMC55"
rotation.dat$PO[rotation.dat$Port.Outside == "SMC65"] <- "SMC55"
rotation.dat$PO[rotation.dat$Port.Outside == "DMC75"] <- "DMC75"

## Port inside
rotation.dat$PI <- NA
rotation.dat$PI[rotation.dat$Port.Inside == "SMC45"] <- "SMC45"
rotation.dat$PI[rotation.dat$Port.Inside == "SMC55"] <- "SMC55"
rotation.dat$PI[rotation.dat$Port.Inside == "SMC65"] <- "SMC55"
rotation.dat$PI[rotation.dat$Port.Inside == "DMC75"] <- "DMC75"

## get a unique net configuration variable
## unique identifier for each set of rotations
rotation.dat$netconfig <- with(rotation.dat, paste(SO, SI, PO, PI, sep = ":"))

rotation.dat$netconfig

##rotation.dat$HAUL <- paste("H", rotation.dat$Haul.., sep = "")
rotation.dat$HAUL <- paste("H", rotation.dat$Haul, sep = "")
rotation.dat$fHAUL <- factor(rotation.dat$HAUL, levels = unique(rotation.dat$HAUL))


## ----eval = TRUE---------------------------------------------------------
## get count per length bin per haul by mesh size
## using the reshape package (makes it easier to process data)
## This should be the same as a pivot (Haul, Carapace length, COMPARTMENT) in excel
library(reshape)

## variables to keep
vars2keep <- c("COMPARTMENT", "Carapace.length", "fHAUL", "COUNT")

## melt the data frame
smc.melt <- melt(smc.dat[, vars2keep],
                  id = c("COMPARTMENT", "Carapace.length", "fHAUL"))

## re-form the dataframe in required format
smc.cast <- cast(smc.melt, Carapace.length + fHAUL  ~ COMPARTMENT  + variable)
smc.cast <- smc.cast[order(smc.cast$fHAUL, smc.cast$Carapace.length), ]

##
## write.csv(subset(smc.cast, fHAUL == "H1"), file = "../data/haul1.csv", row.names = FALSE)

smc.cast[is.na(smc.cast)] <- 0

## merge in the net position
smc.cast <- merge(smc.cast, rotation.dat[, c("fHAUL", "netconfig")])

## merge in the bulk weights
bulk.weight.melt <- melt(unique(smc.dat[ , c("fHAUL", "COMPARTMENT", "Bulk.weight")]),
                         id = c("COMPARTMENT", "fHAUL"))

bulk.weight.cast <- cast(bulk.weight.melt, fHAUL  ~ COMPARTMENT  + variable)

smc.cast <- merge(smc.cast, bulk.weight.cast)

## show the first few rows
head(smc.cast, 2)

## format the subsampling ratio similarly

## double-check that there are unique raising factors per haul
rf.count <- with(smc.dat, table(fHAUL, SUBSRATIO, COMPARTMENT))
apply(rf.count, 1, FUN = function(x){sum(x>0)}) ## yes

## could also check with
unique(smc.dat[, c("fHAUL", "COMPARTMENT", "SUBSRATIO")])

## convert to sub-sampling ratio as in Celtic Warrior
names(smc.dat)[names(smc.dat) == "SUBSRATIO"] <- "SUBSRATIO"
vars2keep <- c("COMPARTMENT", "fHAUL", "SUBSRATIO")
subs.melt <- melt(unique(smc.dat[, vars2keep]), id = c("COMPARTMENT", "fHAUL"))
subs.cast <- cast(subs.melt, fHAUL  ~ COMPARTMENT + variable)

## merge counts and subsampling ratio back together
smc.cast <- merge(smc.cast, subs.cast, by = "fHAUL", all.x = TRUE)

## show first few lines
## Note that this is how the data look just prior to analysis
head(smc.cast, 2)

## ------------------------------------------------------------------------
## Extract the matrix of counts
## Essentialy
count.vars <- c("SMC45_COUNT", "SMC55_COUNT", "SMC65_COUNT", "DMC75_COUNT")

neph.count.mat <- as.matrix(smc.cast[, count.vars])

## Extract the matrix of subsampling ratios
subsratio.vars <- c("SMC45_SUBSRATIO", "SMC55_SUBSRATIO", "SMC65_SUBSRATIO", "DMC75_SUBSRATIO")

subsratio.mat <- as.matrix(smc.cast[, subsratio.vars])

## Create the offset
offset.mat <- log(subsratio.mat / subsratio.mat[,1])

## ----eval = TRUE, fig.allign = "center", fig.cap = "Proportion of Nephrops catch retained per haul. Each point represents the proportion of raised Nephrops catch per haul and length class retained in a given cod-end. The size of the point is proportional to the log of the count.", fig.width = 8, fig.height = 8----

library(ggplot2)

## N.B. to plot the proportions correctly, use the raised counts
## raised count per compartment and haul and carapace length
raised.count.compartment <- neph.count.mat / subsratio.mat

## Get the proportions
prop.compartment <- prop.table(raised.count.compartment, margin = 1)

m <- dim(prop.compartment)[1]

## make a dataframe of the proportions for ggplot
prop.compartment.df <- data.frame(
                         proportion = c(prop.compartment),
                         count = c(raised.count.compartment),
                         Carapace.length = rep(smc.cast$Carapace.length, times = 4),
                         COMPARTMENT = rep(c("SMC45", "SMC55",
                           "SMC65", "DMC75"), each = m),
						fHAUL = rep(smc.cast$fHAUL, times = 4)
                         )

theme_set(theme_bw())

p <- ggplot(prop.compartment.df, aes(x = Carapace.length, y = proportion)) +
  geom_point(col = "red", alpha = 0.2, aes(size = 2*log(count))) +
  facet_wrap(~ COMPARTMENT) + ylab("Proportion of Nephrops per compartment") +
  theme(legend.position = "bottom")

p

p1 <- ggplot(prop.compartment.df, aes(x = Carapace.length, y = proportion)) +
  ##geom_point(colour = "#F8766D", alpha = 0.2, aes(size = 2*log(count))) +
  geom_point(colour = "red", alpha = 0.2, aes(size = 2*log(count))) +
  ##facet_wrap(~ COMPARTMENT) + ylab("Proportion of Nephrops per compartment") +
  facet_grid(fHAUL ~ COMPARTMENT) + ylab("Proportion of Nephrops per compartment") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75)) +
  theme(legend.position = "bottom")

png("../tex/figures/our_lass_by_haul_proportions.png", height = 10, width = 7, res = 300, units = "in")
p1
dev.off()

## ----eval = TRUE, warning = FALSE----------------------------------------
library(nnet)

## First fit is constant proportions
## not accounting for length
## overall proportion
mnom0 <- multinom(neph.count.mat ~ 1 + offset(offset.mat))

## second fit include net configuration/rotations
mnom0.1 <- multinom(neph.count.mat ~ netconfig +
                    offset(offset.mat), data = smc.cast)

## get some predictions from this model
Xpred <- rbind(c(1, 0, 0, 0),
			  c(1, 1, 0, 0),
               c(1, 0, 1, 0),
               c(1, 0, 0, 1),
               c(1, 1/3, 1/3, 1/3))

beta <- cbind(0, t(coef(mnom0.1)))
eta <- Xpred %*% beta
pred.p <- exp(eta) / rowSums(exp(eta))
rownames(pred.p) <- c(unique(smc.cast$netconfig), "Overall")
colnames(pred.p) <- c("CTRL","NSG1","NSG2","SG")

## include carapace length polynomials of different complexity
mnom1 <- update(mnom0.1, . ~ . + poly(Carapace.length, 1),
                data = smc.cast)

mnom2 <- update(mnom0.1, . ~ . + poly(Carapace.length, 2),
                data = smc.cast)

mnom3 <- update(mnom0.1, . ~ . + poly(Carapace.length, 3), data = smc.cast)

AIC(mnom0, mnom0.1, mnom1, mnom2, mnom3)
## looks like a quadratic carapace length effect fits best


## ----eval = TRUE, warnings = FALSE---------------------------------------

## get predictions manually
## CIs not defined in multinomial context but let's try

best.model <- mnom2

## fit coefficients
beta.mu <- c(t(coef(best.model)))

## fit coefficient variance covariance matrix
Sigma <- vcov(best.model)

## number of lengths to predict for
nlength <- 50

pred.length <- seq(min(smc.cast$Carapace.length),
                   max(smc.cast$Carapace.length), length = nlength)

## get the polynomial of lengths
polyfun <- poly(smc.cast$Carapace.length, 2)

## model matrix
Xpred <- cbind(1, 1/3, 1/3, 1/3, predict(polyfun, pred.length))

## number of times to resample predictions to get CIs
nresamp <- 1e3
pred.array <- array(NA, dim = c(nlength, 4, nresamp))

## package to draw from multivariate normal
library(mvtnorm)

for(i in 1:nresamp){
  ##print(i)
  beta0 <- matrix(rmvnorm(1, mean = beta.mu, sigma = Sigma),
                  nrow = 3, byrow = TRUE)
  beta <- cbind(0, t(beta0))
  eta <- Xpred %*% beta
  pred.p <- exp(eta) / rowSums(exp(eta))
  pred.array[ , , i] <- pred.p
  rm(pred.p)
}

## mean across samples
pred.mu <- apply(pred.array, c(1, 2), mean)

## upper across samples
pred.upper <- apply(pred.array, c(1, 2), quantile, p = 0.975)

## lower across samples
pred.lower <- apply(pred.array, c(1, 2), quantile, p = 0.025)

## bring all together in a data frame for ggplot
m <- dim(pred.mu)[1]

   
pred.ci.df <- data.frame(
                COMPARTMENT = rep(c("SMC45", "SMC55",
                  "SMC65", "DMC75"), each = m),
                Carapace.length = rep(pred.length, times = 4),
                proportion = c(pred.mu),
                lower = c(pred.lower),
                upper = c(pred.upper))

## ----eval = TRUE, fig.allign = "center", fig.cap = "Proportion of Nephrops catch retained per haul with fitted multinomial model (without a weight effect or random effects) and associated re-sampled intervals. Null hypothesis of equal retention is displayed as the dashed line at 0.25.", fig.width = 8, fig.height = 8----

p + geom_ribbon(data = pred.ci.df, aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = "blue") +
  geom_line(data = pred.ci.df, aes(x = Carapace.length, y = proportion),
            col = "navy", size = 0.5) +
  geom_hline(aes(yintercept = 0.25), linetype = "dashed")


## ------------------------------------------------------------------------

Yobs <- neph.count.mat
n <- nrow(Yobs)
m <- ncol(Yobs)
X <- model.matrix(best.model)
## here is where the bulk weights enter
Xcond <- smc.cast[, c("CTRL_Bulk.weight", "NSG1_Bulk.weight",
                            "NSG2_Bulk.weight", "SG_Bulk.weight")]
## hauls
gps <- as.numeric(smc.cast$fHAUL)
ngp <- length(unique(gps))

## predicted values of weight
## get median bulk weight per net

median.weights <- with(bulk.weight.melt, tapply(value, COMPARTMENT, median))
median.weights
## construct the "conditional" prediction matrix
npred <- nrow(Xpred)
Xcondpred <- matrix(median.weights, nrow = 1)[rep(1, npred), ]

## Output the data to ADMB
## write the data out
datfile <- "../admbre/multinomialme.dat"
cat("# number of observations n \n", n, "\n", file = datfile)
cat("# number of categories m \n", m, "\n", file = datfile, append = TRUE)
cat("# dimension of parameter vector p \n", ncol(X), "\n",
    file = datfile, append = TRUE)
cat("# dimension of conditional variables q \n", ncol(Xcond), "\n",
    file = datfile, append = TRUE)
cat("# number of groups ngp \n", ngp, "\n", file = datfile, append = TRUE)
cat("# response counts Y \n", file = datfile, append = TRUE)
write.table(Yobs, file = datfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)
cat("# model/design matrix X \n", file = datfile, append = TRUE)
write.table(X, file = datfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)
cat("# conditional matrix Xcond \n", file = datfile, append = TRUE)
write.table(Xcond, file = datfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)
cat("# groups \n", gps, "\n", file = datfile, append = TRUE)
cat("# Offset (log(q[i]/q[1])) \n", file = datfile, append = TRUE)
write.table(offset.mat, file = datfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)

## predictions
cat("# number of prediction rows npred \n", npred, "\n",
    file = datfile, append = TRUE)
cat("# prediction matrix Xpred \n", file = datfile, append = TRUE)
write.table(Xpred, file = datfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)
cat("# conditional prediction matrix Xcondpred \n",
    file = datfile, append = TRUE)
write.table(Xcondpred, file = datfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)

## starting values pinfile
beta0.start <- matrix(0, nrow = ncol(X), ncol = ncol(Yobs) - 1)
u0.start <- matrix(0, nrow = ngp, ncol = m - 1)
nchol <- (m-1)*(m)/2
L.start <- diag(m-1)
L.start[upper.tri(L.start)] <- NA
a.start <- as.numeric(na.omit(c(t(L.start))))
##a.start <- rep(0.1, 3)

pinfile <- "../admbre/multinomialme.pin"
cat("# beta0 \n", file = pinfile)
write.table(beta0.start, file = pinfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)
cat("# betacond \n", 0, "\n", file = pinfile, append = TRUE)
##set.seed(10)
cat("# a \n", a.start, "\n", file = pinfile, append = TRUE)
##cat("# a \n", rep(0.1, m-1), "\n", file = pinfile, append = TRUE)
cat("# u0 \n", file = pinfile, append = TRUE)
write.table(u0.start, file = pinfile, append = TRUE,
            col.names = FALSE, row.names = FALSE)

## NOTE FROM HERE THE CODE NEEDS TO BE RUN IN ADMB
## RUN THE MODEL IN ADMB

## ----eval = TRUE, fig.allign = "center", fig.cap = "Proportion of Nephrops catch retained per haul with fitted ADMB multinomial model (with bulk weights set to the mean bulk per compartment). No random effects were estimated here due to convergence issues which will be investigated.", fig.width = 8, fig.height = 8----

coef.admb <- read.table("../admbre/multinomialme.std", header = TRUE)

##
p0 <- dim(X)[2]
beta.hat <- cbind(0, matrix(subset(coef.admb, name == "beta0")$value,
                            nrow = p0, byrow = TRUE))

betacond.hat <- subset(coef.admb, name == "betacond")$value

## predictions from admb
eta.hat.admb <- matrix(subset(coef.admb, name == "etapred")$value,
                       ncol = 4, byrow = TRUE)
P.hat.admb <- exp(eta.hat.admb) / rowSums(exp(eta.hat.admb))
eta.se.admb <- matrix(subset(coef.admb, name == "etapred")$std.dev,
                      ncol = 4, byrow = TRUE)

eta.lwr.admb <- eta.hat.admb + qnorm(0.025) * eta.se.admb
eta.upr.admb <- eta.hat.admb + qnorm(0.975) * eta.se.admb

P.admb.lwr <- exp(eta.lwr.admb) / rowSums(exp(eta.lwr.admb))
P.admb.upr <- exp(eta.upr.admb) / rowSums(exp(eta.upr.admb))

plot.pred.df.admb <- data.frame(
                       COMPARTMENT = rep(c("Control", "Nephrops sorting grid 1",
                           "Nephrops sorting grid 2", "Swedish grid"),
                         each = nrow(P.hat.admb)),
                       Carapace.length = rep(pred.length, times = 4),
                       proportion = c(P.hat.admb),
                       lwr = c(P.admb.lwr),
                       upr = c(P.admb.upr))

p + geom_ribbon(data = plot.pred.df.admb, aes(ymin = lwr, ymax = upr),
                alpha = 0.3, fill = "blue") +
  geom_line(data = plot.pred.df.admb, aes(x = Carapace.length, y = proportion),
            col = "navy", size = 0.5) +
  geom_hline(aes(yintercept = 0.25), linetype = "dashed")


## ----eval = FALSE--------------------------------------------------------
##
## uhat <- matrix(subset(coef.admb, name == "u0")$value, ncol = m-1, byrow = TRUE)
##
## par(mar = c(4, 4, 2, 2))
## matplot(1:12, uhat, pch = c(2, 8, 19), col = 1, xlab = "", ylab = "", xaxt = "n")
## axis(side = 1, at = 1:12, labels = c(1:3, 5:13))
## abline(h = 0, lty = 2)
## mtext(side = 1, line = 2.5, text = "Haul number")
## mtext(side = 2, line = 2.5, text = "Log-odds random effects")
## legend("topright", legend = c("NSG1/CTRL", "NSG2/CTRL", "SG/CTRL"), pch = c(2, 8, 19))
##
## plot(as.data.frame(uhat)) ## not as strong correlation int he random effects
##
## ## Correlation of the random effects
## a.hat <- subset(coef.admb, name == "a")$value
##
## ii <- 1
##
## L <- matrix(0, 3, 3)
## for(i in 1:3){
##   for(j in 1:i){
##     L[i,j] <- a.hat[ii]
##     ii <- ii + 1
##   }
## }
##
## sigma.hat <- L %*% t(L)
##
## cov2cor(sigma.hat)
## cor(data.frame(uhat))
##

## ----eval = FALSE--------------------------------------------------------
## ## Haul level predictions (nest steps)
##

