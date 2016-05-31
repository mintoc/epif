##------------------------------------
## Compare betabin with glmm
## CM, BB: Thu Jul 23 2015
## NOTE: Sort out deviance calculations mathematically
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

sub.dat <- subset(neph.dat, Mesh.Size %in% c("70mm", "100mm"))

sub.dat <- reshape(sub.dat,
                timevar = "Mesh.Size",
                idvar = c("fHAUL", "Carapace.Length"),
                direction = "wide",
                drop = names(sub.dat)[!names(sub.dat) %in% c("Mesh.Size", "fHAUL", "Carapace.Length", "fMesh.Size", "SUBSRATIO", "COUNT")])

sub.dat$COUNT.70mm[is.na(sub.dat$COUNT.70mm)] <- 0
sub.dat$COUNT.100mm[is.na(sub.dat$COUNT.100mm)] <- 0


## fill in the SUBSRATIO
subs.uniq <- unique(neph.dat[, c("SUBSRATIO","fHAUL", "Mesh.Size")])

sub.dat$SUBSRATIO.70mm <- merge(sub.dat, subset(subs.uniq, Mesh.Size == "70mm"), by = c("fHAUL"), sort = FALSE)$SUBSRATIO

sub.dat$SUBSRATIO.100mm <- merge(sub.dat, subset(subs.uniq, Mesh.Size == "100mm"), by = c("fHAUL"), sort = FALSE)$SUBSRATIO

## Have a look at the data
with(sub.dat, plot(Carapace.Length, (COUNT.70mm/SUBSRATIO.70mm)/(COUNT.70mm/SUBSRATIO.70mm + COUNT.100mm/SUBSRATIO.100mm)))

## DROP REALLY BIG FOR PRELIM ANALYSIS
sub.dat <- subset(sub.dat, Carapace.Length < 60)

sub.dat$TOTAL <- with(sub.dat, COUNT.70mm/SUBSRATIO.70mm + COUNT.100mm/SUBSRATIO.100mm)

with(sub.dat, plot(Carapace.Length, (COUNT.70mm/SUBSRATIO.70mm)/(COUNT.70mm/SUBSRATIO.70mm + COUNT.100mm/SUBSRATIO.100mm), cex = 1/2 * log(TOTAL), pch = 19, col = "#FF000020"))

## SOME MODELS
## BASIC BINOMIAL
## NOTE: Confirm offset mathematically not at 5.30
library(mgcv)

binom.70.100 <- glm(cbind(COUNT.70mm, COUNT.100mm) ~ poly(Carapace.Length, 3) + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)), data = sub.dat, family = binomial)

binom.70.100.gam <- gam(cbind(COUNT.70mm, COUNT.100mm) ~ s(Carapace.Length) + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)), data = sub.dat, family = binomial)

curve(dchisq(x, df = binom.70.100$df.residual), from = 0, to = 2e3, n = 1e3)

## OVERLAY FIT
with(sub.dat, plot(Carapace.Length, (COUNT.70mm/SUBSRATIO.70mm)/(COUNT.70mm/SUBSRATIO.70mm + COUNT.100mm/SUBSRATIO.100mm), cex = 1/2 * log(TOTAL), pch = 19, col = "#FF000020"))

curve(predict(binom.70.100, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "blue")

curve(predict(binom.70.100.gam, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "blue", lty = 2)

##------
## GLMM
##------
library(lme4)
binom.70.100.glmm <- glmer(cbind(COUNT.70mm, COUNT.100mm) ~ poly(Carapace.Length, 3) + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)) + (1 | fHAUL), data = sub.dat, family = binomial)

##------
## GAMM 
##------
binom.70.100.gamm <- gamm(cbind(COUNT.70mm, COUNT.100mm) ~ s(Carapace.Length), random = list(fHAUL = ~ 1),  data = sub.dat, family = binomial)

## check but looks like deviance is
sum(resid(binom.70.100.gamm$lme)^2)
## so still high 

## what do these look like?
## OVERLAY FIT
with(sub.dat, plot(Carapace.Length, (COUNT.70mm/SUBSRATIO.70mm)/(COUNT.70mm/SUBSRATIO.70mm + COUNT.100mm/SUBSRATIO.100mm), cex = 1/2 * log(TOTAL), pch = 19, col = "#FF000020"))

## glm
curve(predict(binom.70.100, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "blue")

## gam
curve(predict(binom.70.100.gam, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "blue", lty = 2)

## glmm
curve(predict(binom.70.100.glmm, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response", re.form = NA), add = TRUE, lwd = 2, col = "purple")

curve(predict(binom.70.100.gamm$gam, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "purple", lty = 2)

##---------------
## BETA-BINOMIAL 
##---------------
library(gamlss)

vars2keep <- c("COUNT.70mm", "COUNT.100mm", "Carapace.Length", "SUBSRATIO.70mm", "SUBSRATIO.100mm")

binom.70.100.gamlss <- gamlss(cbind(COUNT.70mm, COUNT.100mm) ~ poly(Carapace.Length, 3) + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)), family = BI, data = sub.dat[, vars2keep])

betabinom.70.100.gamlss <- gamlss(cbind(COUNT.70mm, COUNT.100mm) ~ poly(Carapace.Length, 3) + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)), family = BB, data = sub.dat[, vars2keep])

##-----
## GEE 
##-----
library(geepack)
binom.70.100.gee <- geeglm(cbind(COUNT.70mm, COUNT.100mm) ~ poly(Carapace.Length, 3) + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)), id = fHAUL, family = binomial, data = sub.dat, corstr = "exchangeable")

## OVERLAY FIT
with(sub.dat, plot(Carapace.Length, (COUNT.70mm/SUBSRATIO.70mm)/(COUNT.70mm/SUBSRATIO.70mm + COUNT.100mm/SUBSRATIO.100mm), cex = 1/2 * log(TOTAL), pch = 19, col = "#FF000020"))

## glmm
curve(predict(binom.70.100.glmm, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response", re.form = NA), add = TRUE, lwd = 2, col = "purple")

## betabin
curve(predict(betabinom.70.100.gamlss, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "purple", lty = 2)

## gee
curve(predict(binom.70.100.gee, newdata = data.frame(Carapace.Length = x, SUBSRATIO.70mm = 1, SUBSRATIO.100mm = 1), type = "response"), add = TRUE, lwd = 2, col = "purple", lty = 4)

##----------------
## COEF ESTIMATES 
##----------------
## plot the fixed effect estimates

coef.names <- c("intercept", paste("CL.poly", 1:3, sep = ""))
glm.par <- as.data.frame(summary(binom.70.100)$coefficients)
glm.par$coef <- coef.names
glm.par$method <- "GLM"
glm.par$coef.num <- 1:4
##
glmm.par <- as.data.frame(summary(binom.70.100.glmm)$coefficients)
glmm.par$coef <- coef.names
glmm.par$method <- "GLMM"
glmm.par$coef.num <- 1:4
##
betabin.par <- summary(betabinom.70.100.gamlss)
betabin.par <- as.data.frame(betabin.par[-5, ])
betabin.par$coef <- coef.names
betabin.par$method <- "BB"
betabin.par$coef.num <- 1:4
##
gee.par <- summary(binom.70.100.gee)$coefficients
gee.par$coef <- coef.names
gee.par$method <- "GEE"
names(gee.par)[2] <- "Std. Error"
gee.par$coef.num <- 1:4

all.par <- rbind(glm.par[, c(1:2, 5:7)],
                 glmm.par[, c(1:2, 5:7)],
                 betabin.par[, c(1:2, 5:7)],
                 gee.par[, c(1:2, 5:7)])

all.par$upr <- all.par$Estimate + 2 * all.par[, "Std. Error"]
all.par$lwr <- all.par$Estimate - 2 * all.par[, "Std. Error"]

library(ggplot2)
pd <- position_dodge(width=0.2,height=NULL)

pdf("../tex/figures/binomial_coef_comparison.pdf", height = 6, width = 7)
p <- ggplot(all.par, aes(coef.num, Estimate, color = method)) +
  geom_point(aes(shape=method),size=4, position=pd) +
  scale_shape_manual(name = "Method",values = c(15, 16, 17, 18)) +
  scale_color_manual(name="Method",values=c("coral","steelblue", "turquoise", "magenta")) +
  theme_bw() +
  scale_x_continuous("Coefficient", breaks=1:4, labels = glm.par$coef) +
  ##scale_y_continuous(trans = 'asinh') + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.1, position=pd)
print(p)
dev.off()

##---------
## SANDBOX
##---------
## NEED TO FIGURE OUT DEVIANCE CALCS FROM FIRST PRINCIPLES
binom.70.100.null <- glm(cbind(COUNT.70mm, COUNT.100mm) ~ 1 + offset(log(SUBSRATIO.70mm/SUBSRATIO.100mm)), data = sub.dat, family = binomial)

null.like <- with(sub.dat, dbinom(COUNT.70mm, size = COUNT.70mm + COUNT.100mm, p = predict(binom.70.100.null, type = "response")))

alt.like <- with(sub.dat, dbinom(COUNT.70mm, size = COUNT.70mm + COUNT.100mm, p = predict(binom.70.100, type = "response")))
