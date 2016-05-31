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

## load the model fits
load("best_fit.RData") ## opt and rep

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

prop.plot <-
  ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion, colour = fHAUL)) +
  geom_point(aes(size = pch.cex)) + scale_size_continuous(range = c(1,4.5)) +
  facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75)) + theme(legend.position = "none")

theme_set(theme_bw(base_size = 10))

png("../tex/ICES paper review/figures/Figure2.png", height = 15, width = 17, units = "cm", res = 500)
prop.plot
dev.off()

## Plot to add overall fit back in
prop.plot.nocolour <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion)) +
  ##geom_point(aes(size = pch.cex), colour = "#F8766D") +
  geom_point(aes(size = pch.cex), pch = 1, colour = "darkgrey") +
  scale_size_continuous(range = c(1,4.5)) +
facet_wrap(~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end") + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75)) + theme(legend.position = "none")

## Plot to add by-haul fits in
## By Haul
prop.plot.haul <- ggplot(prop.mesh.df, aes(x = Carapace.length, y = proportion)) +
  geom_point(aes(size = pch.cex), shape = 1, colour = "darkgrey") +
  scale_size_continuous(range = c(1,4.5)) +
  facet_grid(fHAUL ~ Mesh.Size) + ylab("Proportion of raised Nephrops per cod-end")  + xlab("Carapace length (mm)") + scale_y_continuous(breaks=c(0.25, 0.5, 0.75))
##scale_colour_gradientn(colours=rainbow(4))

##---------------
## INCLUDE FITS
##---------------
## read back in the random effects results
m <- 4
coef.tmb <- summary(rep)
uhat <- matrix(coef.tmb[rownames(coef.tmb) == "u0", 1], ncol = m-1)

png("../tex/ICES paper review/figures/Figure3.png", height = 15, width = 17, units = "cm", res = 500)
par(mar = c(4, 4, 2, 2))
matplot(1:12, uhat, pch = c(2, 8, 19), col = 1, xlab = "", ylab = "", xaxt = "n")
axis(side = 1, at = 1:12)
abline(h = 0, lty = 2)
mtext(side = 1, line = 2.5, text = "Haul number")
mtext(side = 2, line = 2.5, text = "Log-odds random effects")
legend("topright", legend = c("80mm/70mm", "90mm/70mm", "100mm/70mm"), pch = c(2, 8, 19))
dev.off()

## Correlation of the random effects
a.hat <- coef.tmb[rownames(coef.tmb) == "a", 1]

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

write.csv(x = round(sigma.hat, 2), "../data/sigma_u_hat.csv")
write.csv(x = round(cov2cor(sigma.hat), 2), "../data/correlation_u_hat.csv")

## PREDICTIONS
npred <- 50
CLpred <- seq(min(our.lass.neph.cast$Carapace.length), max(our.lass.neph.cast$Carapace.length), length = npred)

eta.hat.tmb <- matrix(coef.tmb[rownames(coef.tmb) == "etapred", 1], ncol = 4)
eta.se.tmb <- matrix(coef.tmb[rownames(coef.tmb) == "etapred", 2], ncol = 4)

eta.lwr.tmb <- eta.hat.tmb + qnorm(0.025) * eta.se.tmb
eta.upr.tmb <- eta.hat.tmb + qnorm(0.975) * eta.se.tmb

P.hat.tmb <- exp(eta.hat.tmb) / rowSums(exp(eta.hat.tmb))
P.tmb.lwr <- exp(eta.lwr.tmb) / rowSums(exp(eta.lwr.tmb))
P.tmb.upr <- exp(eta.upr.tmb) / rowSums(exp(eta.upr.tmb))

plot.pred.df.tmb <- data.frame(
                       Mesh.Size = factor(rep(colnames(Y), each = length(CLpred))),
                       Carapace.length = rep(CLpred, times = 4),
                       proportion = c(P.hat.tmb),
                       lwr = c(P.tmb.lwr),
                       upr = c(P.tmb.upr))

levels.switch <- c("m70mm_Count" = "70mm", "m80mm_Count" = "80mm", "m90mm_Count" = "90mm", "m100mm_Count" = "100mm")

levels(plot.pred.df.tmb$Mesh.Size) <- levels.switch[levels(plot.pred.df.tmb$Mesh.Size)]

png("../tex/ICES paper review/figures/Figure5.png", height = 15, width = 17, units = "cm", res = 500)
prop.plot.nocolour +
  ##geom_ribbon(data = plot.pred.df.tmb, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.4) +
  geom_line(data = plot.pred.df.tmb, aes(y = proportion), linetype = "solid") +
  geom_line(data = plot.pred.df.tmb, aes(y = lwr), linetype = "dashed") +
  geom_line(data = plot.pred.df.tmb, aes(y = upr), linetype = "dashed")
  ##geom_line(data = plot.pred.df.tmb, aes(y = dm.proportion), colour = "blue")
dev.off()

## HAUL LEVEL PREDICTIONS
## predictions at the haul level with and without random effects
## note haul levels ordered by group
## so reffs ordered
levels(our.lass.neph.cast$fHAUL)

par.hat <- rep$par.fixed
fixed.names <- names(par.hat)

## beta
##beta0.hat <- matrix(par.hat[fixed.names == "beta0"], nr = 2, nc = 3, byrow = FALSE)
beta0.hat <- matrix(0, nr = 2, nc = 3, byrow = FALSE)
## use the mapping
beta0.hat[c(1,2,3,5)] <- par.hat[fixed.names == "beta0"]
##beta0.hat[1,] <- par.hat[fixed.names == "beta0"][c(1,3,4)]
##beta0.hat[2,1] <- par.hat[fixed.names == "beta0"][2]
beta.hat <- cbind(0, beta0.hat)
CL <- our.lass.neph.cast$Carapace.length
CL.sc <- with(our.lass.neph.cast, (CL - mean(Carapace.length)) / sd(Carapace.length))
X <- cbind(1, CL.sc)

## NET POSITION
betacond.hat <- par.hat[fixed.names == "betacond"]
##NP <- betacond.hat[1] * PI + betacond.hat[2] * SI + betacond.hat[3] * SO
NP <- betacond.hat[1] * PI  + betacond.hat[2] * SO

## WEIGHTS
gammaw.hat <- diag(par.hat[fixed.names == "gammaw"])

bulk.weights.mat <- as.matrix(our.lass.neph.cast[, c("m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")])
colnames(bulk.weights.mat) <- rownames(bulk.weights.mat) <- NULL
W <- scale.fun(bulk.weights.mat) ## to maintain relative weights of nets

## WEIGHTS x CL INTERACTION
WxCL <- scale.fun(bulk.weights.mat * our.lass.neph.cast[, "Carapace.length"])

## NET POSITION X CL
##NPxCL <-
##  betacond.hat[4] * PI * CL.sc +
##  betacond.hat[5] * SI * CL.sc +
##  betacond.hat[6] * SO * CL.sc

NPxCL <-
  betacond.hat[3] * PI * CL.sc +
  betacond.hat[4] * SO * CL.sc


## WEIGHT x CL
WxCL <- scale.fun(bulk.weights.mat * our.lass.neph.cast[, "Carapace.length"])

##gammawcl.hat <- diag(par.hat[fixed.names == "gammawcl"])
gammawcl.hat <- diag(c(par.hat[fixed.names == "gammawcl"], 0, 0))

eta.hat <- X %*% beta.hat + NP + W %*% gammaw.hat + NPxCL + WxCL %*% gammawcl.hat + cbind(0, uhat)[gps,]

P.hat <- exp(eta.hat) / rowSums(exp(eta.hat))

## output parameter table from random effects fit
## Note won't work with non-significant parameters removed
## output this with zeroes in place of removed parameters
par.tab <- data.frame(Parameter = c("Intercept", "Carapace length (CL)", "Port inside (PI)", "Starboard inside (SI)", "Starboard outside (SO)", "PI x CL", "SI x CL", "SO x CL", "Cod-end weight (W)*", "W x CL*"), "70/70" = 0, "80/70" = 0, "90/70" = 0, "100/70" = 0)
names(par.tab) <- c("Parameter", "70/70", "80/70", "90/70", "100/70")
par.tab[1:2, 3:5] <- paste(apply(round(coef.tmb[rownames(coef.tmb) == "beta0", ], 3), 1, paste, collapse = "("), ")", sep = "")
par.tab[3:8, 2] <- paste(apply(round(coef.tmb[rownames(coef.tmb) == "betacond", ], 3), 1, paste, collapse = "("), ")", sep = "")
par.tab[9, 2:5] <- paste(apply(round(coef.tmb[rownames(coef.tmb) == "gammaw", ], 3), 1, paste, collapse = "("), ")", sep = "")
par.tab[10, 2:5] <- paste(apply(round(coef.tmb[rownames(coef.tmb) == "gammawcl", ], 3), 1, paste, collapse = "("), ")", sep = "")

write.csv(par.tab, file = "../data/parameter_table.csv", row.names = FALSE)

## predictions at haul level for model with no random effects
load("best_fit_nore.RData")

par.hat <- rep_nore$par.fixed
fixed.names <- names(par.hat)

## beta
##beta0.hat <- matrix(par.hat[fixed.names == "beta0"], nr = 2, nc = 3, byrow = FALSE)
beta0.hat <- matrix(0, nr = 2, nc = 3, byrow = FALSE)
## use the mapping
beta0.hat[c(1,2,3,5)] <- par.hat[fixed.names == "beta0"]
##beta0.hat[1,] <- par.hat[fixed.names == "beta0"][c(1,3,4)]
##beta0.hat[2,1] <- par.hat[fixed.names == "beta0"][2]
beta.hat <- cbind(0, beta0.hat)
CL <- our.lass.neph.cast$Carapace.length
CL.sc <- with(our.lass.neph.cast, (CL - mean(Carapace.length)) / sd(Carapace.length))
X <- cbind(1, CL.sc)

## NET POSITION
betacond.hat <- par.hat[fixed.names == "betacond"]
##NP <- betacond.hat[1] * PI + betacond.hat[2] * SI + betacond.hat[3] * SO
NP <- betacond.hat[1] * PI  + betacond.hat[2] * SO

## WEIGHTS
gammaw.hat <- diag(par.hat[fixed.names == "gammaw"])

bulk.weights.mat <- as.matrix(our.lass.neph.cast[, c("m70mm_Total.catch", "m80mm_Total.catch", "m90mm_Total.catch", "m100mm_Total.catch")])
colnames(bulk.weights.mat) <- rownames(bulk.weights.mat) <- NULL
W <- scale.fun(bulk.weights.mat) ## to maintain relative weights of nets

## WEIGHTS x CL INTERACTION
WxCL <- scale.fun(bulk.weights.mat * our.lass.neph.cast[, "Carapace.length"])

## NET POSITION X CL
##NPxCL <-
##  betacond.hat[4] * PI * CL.sc +
##  betacond.hat[5] * SI * CL.sc +
##  betacond.hat[6] * SO * CL.sc

NPxCL <-
  betacond.hat[3] * PI * CL.sc +
  betacond.hat[4] * SO * CL.sc


## WEIGHT x CL
WxCL <- scale.fun(bulk.weights.mat * our.lass.neph.cast[, "Carapace.length"])

##gammawcl.hat <- diag(par.hat[fixed.names == "gammawcl"])
gammawcl.hat <- diag(c(par.hat[fixed.names == "gammawcl"], 0, 0))

eta.hat.nore <- X %*% beta.hat + NP + W %*% gammaw.hat + NPxCL + WxCL %*% gammawcl.hat

P.hat.nore <- exp(eta.hat.nore) / rowSums(exp(eta.hat.nore))

pred.df <- data.frame(proportion = as.numeric(unlist(P.hat)),
                      proportion.nore = as.numeric(unlist(P.hat.nore)),
                      Carapace.length = CL,
                      fHAUL = paste("H", gps, sep = ""),
                      Mesh.Size = rep(count.vars, each =  nrow(X)))

levels(pred.df$Mesh.Size) <- levels.switch[levels(pred.df$Mesh.Size)]

png("../tex/ICES paper review/figures/Figure4.png", height = 25, width = 17, units = "cm", res = 500)
prop.plot.haul + geom_line(data = pred.df) + theme(legend.position = "none") + geom_line(data = pred.df, aes(y = proportion.nore), linetype = "dashed", colour = "black")
dev.off()

## odds-ratio predictions
lodds.mat <- coef.tmb[grep("lodds", rownames(coef.tmb)), ]

odds.df <- data.frame(comparison = rownames(lodds.mat),
                      CL = CLpred,
                      lodds = lodds.mat[, "Estimate"],
                      lodds.se = lodds.mat[, "Std. Error"],
                      odds = exp(lodds.mat[, "Estimate"]))

## odds.df$odds2 <- c(P.hat.tmb[,1] / P.hat.tmb[,2],
##                    P.hat.tmb[,1] / P.hat.tmb[,3],
##                    P.hat.tmb[,1] / P.hat.tmb[,4],
##                    P.hat.tmb[,2] / P.hat.tmb[,3],
##                    P.hat.tmb[,2] / P.hat.tmb[,4],
##                    P.hat.tmb[,3] / P.hat.tmb[,4])

odds.df$comparison <- factor(as.character(odds.df$comparison), levels = c("lodds_70_80", "lodds_70_90", "lodds_70_100", "lodds_80_90", "lodds_80_100", "lodds_90_100"))

odds.df$comparison2 <- NA
odds.df$comparison2[odds.df$comparison == "lodds_70_80"] <- "p70 / p80"
odds.df$comparison2[odds.df$comparison == "lodds_70_90"] <- "p70 / p90"
odds.df$comparison2[odds.df$comparison == "lodds_70_100"] <- "p70 / p100"
odds.df$comparison2[odds.df$comparison == "lodds_80_90"] <- "p80 / p90"
odds.df$comparison2[odds.df$comparison == "lodds_80_100"] <- "p80 / p100"
odds.df$comparison2[odds.df$comparison == "lodds_90_100"] <- "p90 / p100"

odds.df$comparison2 <- factor(odds.df$comparison2, levels = c("p70 / p80", "p70 / p90", "p70 / p100", "p80 / p90", "p80 / p100", "p90 / p100"))

alpha <- 0.05 ## signficance level
## regular CI
qnorm(1 - alpha / 2)
## corrected for 6 comparisons
(q.val <- qnorm(1 - alpha / (2 * 6)))
##(q.val <- qnorm(1 - alpha / 2))

odds.df$odds.upr <- with(odds.df, exp(lodds + q.val * lodds.se))
odds.df$odds.lwr <- with(odds.df, exp(lodds - q.val * lodds.se))

png("../tex/ICES paper review/figures/Figure6.png", height = 12, width = 17, units = "cm", res = 500)
ggplot(odds.df, aes(x = CL, y = odds)) + facet_wrap(~ comparison2) + geom_ribbon(aes(ymin = odds.lwr, ymax = odds.upr), fill = "grey") + geom_line() + geom_hline(yintercept = 1, linetype = "dashed") + xlab("Carapace length (mm)") + ylab("Odds ratio")##+ geom_line(aes(y = odds2), col = "blue", linetype = "dashed")
dev.off()

## p80
p80 <- coef.tmb[grep("p_80_70", rownames(coef.tmb)), ]

p80.upr <- p80[, "Estimate"] + 1.96 * p80[, "Std. Error"]
p80.lwr <- p80[, "Estimate"] - 1.96 * p80[, "Std. Error"]

matplot(CLpred, cbind(p80[, "Estimate"], p80.lwr, p80.upr), lty = c(1, 2, 2), type = "l", col = 1, ylim = c(0,1))
abline(h = 0.5, col = "grey")

##------------------
## ALL EFFECTS PLOT HERE!
##------------------

load("best_fit.RData") ## opt and rep

par.hat <- rep$par.fixed
fixed.names <- names(par.hat)

## beta
##beta0.hat <- matrix(par.hat[fixed.names == "beta0"], nr = 2, nc = 3, byrow = FALSE)
beta0.hat <- matrix(0, nr = 2, nc = 3, byrow = FALSE)
## use the mapping
beta0.hat[c(1,2,3,5)] <- par.hat[fixed.names == "beta0"]

##beta0.hat[1,] <- par.hat[fixed.names == "beta0"][c(1,3,4)]
##beta0.hat[2,1] <- par.hat[fixed.names == "beta0"][2]
beta.hat <- cbind(0, beta0.hat)

betacond.hat <- par.hat[fixed.names == "betacond"]

npred <- 10
cl <- seq(min(our.lass.neph.cast$Carapace.length), max(our.lass.neph.cast$Carapace.length), length = npred)

gammaw <- par.hat[fixed.names == "gammaw"]
gammawcl <- par.hat[fixed.names == "gammawcl"]

## weights
W.names <- paste("m", seq(70,100, by = 10), "mm_Total.catch", sep = "")
weights <- unlist(unique(our.lass.neph.cast[, W.names]))
weight.values <- as.numeric(quantile(weights, p = c(0.1, 0.5, 0.9)))

## observed rotations were:
## Hauls 1-3: 70,90,80,100
## Hauls 4-6: 90,70,100,80
## Hauls 7-9: 100,80,70,90
## Hauls 10-12: 80,100,90,70

## use the actual rotations of the trial

rotations <- list(
    "70,90,80,100" = list(
        pi = c(0, 0, 1, 0),
        si = c(0, 1, 0, 0),
        so = c(0, 0, 0, 1)
        ),
    "90,70,100,80" = list(
        pi = c(1, 0, 0, 0),
        si = c( 0, 0, 0, 1),
        so = c(0, 1, 0, 0)
              ),
    "100,80,70,90" = list(
        pi = c(0, 1, 0, 0),
        si = c(1, 0, 0, 0),
        so = c(0, 0, 1, 0)
              ),
    "80,100,90,70" = list(
        pi = c(0, 0, 0, 1),
        si = c(0, 0, 1, 0),
        so = c(1, 0, 0, 0)
        )
    )

w <- mean(unlist(our.lass.neph.cast[, W.names]))

get.prop <- function(beta, betacond, cl, pi, si, so, gammaw, gammawcl, w, rotation.name){
    n <- length(cl)
    ## carapace lengths
    CL.sc <- with(our.lass.neph.cast, (cl - mean(Carapace.length)) / sd(Carapace.length))
    X <- cbind(1, CL.sc)
    ## NET POSITION
    PI <- matrix(pi, nr = 1, nc = 4)[rep(1, n),]
    SI <- matrix(si, nr = 1, nc = 4)[rep(1, n),]
    SO <- matrix(so, nr = 1, nc = 4)[rep(1, n),]
    ##NP <- betacond[1] * PI + betacond[2] * SI + betacond[3] * SO
    NP <- betacond[1] * PI + betacond[2] * SO
    ## WEIGHTS
    gammaw.hat <- diag(gammaw)
    W.names <- paste("m", seq(70,100, by = 10), "mm_Total.catch", sep = "")
    mean.W <- mean(unlist(our.lass.neph.cast[, W.names]))
    sd.W <- sd(unlist(our.lass.neph.cast[, W.names]))
    Wsc.means <- (w - mean.W) / sd.W
    W <- matrix(as.numeric(Wsc.means), nr = 1, nc = 4)[rep(1, n), ]
    ## NET POSITION X CL
    ##NPxCL <-
        ##betacond.hat[4] * PI * CL.sc +
            ##betacond.hat[5] * SI * CL.sc +
                ##betacond.hat[6] * SO * CL.sc
    NPxCL <-
      betacond.hat[3] * PI * CL.sc +
            betacond.hat[4] * SO * CL.sc
    ## WEIGHT x CL
    mean.WxCL <- mean(unlist(our.lass.neph.cast[, W.names] * our.lass.neph.cast[, "Carapace.length"]))
    sd.WxCL <- sd(unlist(our.lass.neph.cast[, W.names] * our.lass.neph.cast[, "Carapace.length"]))
    WxCL <- (t(sapply(cl, FUN = function(x){ matrix(x * matrix(rep(w, 4), nrow = 1))})) - mean.WxCL) / sd.WxCL
    ##gammawcl.hat <- diag(gammawcl)
    gammawcl.hat <- diag(c(gammawcl,0,0))
    ##
    eta <- X %*% beta.hat + NP + W %*% gammaw.hat + NPxCL + WxCL %*% gammawcl.hat
    p <- exp(eta) / rowSums(exp(eta))
    p.df <- data.frame(Carapace.length = rep(cl, times = 4),
                       proportion = c(p),
                       Mesh.Size = rep(c("70mm", "80mm", "90mm", "100mm"), each = n),
                       weight = w,
                       rotation = rotation.name)
    return(p.df)
}

## temporary for interpreting the parameter table
tmp <- NULL
for(i in 1:3){
    tmp0 <- get.prop(beta = beta.hat, betacond = betacond.hat, cl = cl, pi = rep(1/4,4), si = rep(1/4,4), so = rep(1/4,4), gammaw = gammaw, gammawcl = gammawcl, w = weight.values[i], rotation.name = "test")
    tmp <- rbind(tmp, tmp0)
}

tmp$Mesh.Size <- factor(as.character(tmp$Mesh.Size), levels = c("70mm", "80mm", "90mm", "100mm"))
tmp$Weight <- factor(ifelse(tmp$weight == weight.values[1], "Low",
                       ifelse(tmp$weight == weight.values[2], "Medium", "High")), levels = c("Low", "Medium", "High"))

ggplot(tmp, aes(x = Carapace.length, y = proportion)) + facet_wrap( ~ Mesh.Size, ncol = 4) + geom_line(aes(linetype = Weight)) + xlab("Carapace length (mm)") + ylab("Proportion retained") + ylim(0,1) + scale_linetype_manual(values = c("dotted", "dashed", "solid")) + theme(legend.position = "bottom")

p.df <- NULL
for(i in 1:length(rotations)){
    for(j in 1:length(weight.values)){
        p.tmp <- get.prop(beta = beta.hat, betacond = betacond.hat, cl = cl, pi = rotations[[i]]$pi, si = rotations[[i]]$si, so = rotations[[i]]$so, gammaw = gammaw, gammawcl = gammawcl, w = weight.values[j], rotation.name = names(rotations)[i])
        p.df <- rbind(p.df, p.tmp)
    }
}

p.df$Mesh.Size <- factor(as.character(p.df$Mesh.Size), levels = c("70mm", "80mm", "90mm", "100mm"))
p.df$Weight <- factor(ifelse(p.df$weight == weight.values[1], "Low",
                       ifelse(p.df$weight == weight.values[2], "Medium", "High")), levels = c("Low", "Medium", "High"))

png("../tex/ICES paper review/figures/Figure7.png", height = 15, width = 17, units = "cm", res = 500)
ggplot(p.df, aes(x = Carapace.length, y = proportion, group = Weight)) + facet_grid(rotation ~ Mesh.Size) + geom_line(aes(linetype = Weight)) + xlab("Carapace length (mm)") + ylab("Proportion retained") + ylim(0,1) + scale_linetype_manual(values = c("dotted", "dashed", "solid")) + theme(legend.position = "bottom")
dev.off()

##---------
## SANDBOX - effect of weight
##---------



prop.plot.nocolour + geom_line(data = p.df)

matplot(CL, p, type = "l", lty = 1)

## from tmb
eta.hat <- matrix(summary(rep)[rownames(summary(rep)) == "etapred", 1], nc = 4)

p2 <- exp(eta.hat) / rowSums(exp(eta.hat))

matplot(CL, p, type = "l", lty = 1)
matlines(CL, p2, type = "l", lty = 2, lwd = 2)
