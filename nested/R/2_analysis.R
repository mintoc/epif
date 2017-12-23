##-----------------------
## Code to implement nested catch comparison analysis
## CM: August 4 2017
## Notes:
##-----------------------

library(TMB)
library(mgcv)
library(ggplot2); theme_set(theme_bw())
all.cast <- read.csv("../data/dual_data_cast.csv")

## subset for a species
spp.vec <- unique(all.cast$Species)

##for(i in 1:length(spp.vec)){
##i <- 3
print(i)
spp <- spp.vec[i]
dat <- subset(all.cast, Species == spp)
N <- dat[, c("Control_Count", "Lower.codend_Count", "Upper.codend_Count", "Test_Count")]
S <- dat[, c("Control_SUBSRATIO", "Lower.codend_SUBSRATIO", "Upper.codend_SUBSRATIO", "Test_SUBSRATIO")]
## level 1
O1 <- matrix(log(S[, "Test_SUBSRATIO"] / S[, "Control_SUBSRATIO"]))
N1 <- as.matrix(N[, c("Test_Count", "Control_Count")])
O2 <- matrix(log(S[, "Lower.codend_SUBSRATIO"] / S[, "Upper.codend_SUBSRATIO"]))
N2 <- as.matrix(N[, c("Lower.codend_Count", "Upper.codend_Count")])

## set up the splines
x1 <- min(dat$Length)
x2 <- max(dat$Length)

## See Eilers and Marx 2010
tpower <- function(x, t, p){
    ## Truncated p-th power function
    (x - t) ^ p * (x > t)
}

bbase <- function(x, xl, xr, ndx, deg){
    ## Construct a B-spline basis of degree ’deg’
    dx <- (xr - xl) / ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B
}
Z1 <- Z2 <- bbase(x = dat$Length, xl = x1, xr = x2, ndx = 20, deg = 3)

## knots for plotting
xl = x1; xr = x2; ndx = 20; deg = 3
dx <- (xr - xl) / ndx
knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)

Y1 <- N1
Y2 <- N2
X1 <- X2 <- cbind(1, as.numeric(dat$Day.night == "Night")) ##matrix(1, nrow = nrow(Y1))
L1pred <- L2pred <- seq(x1, x2, length = 100)
X1pred <- X2pred <- cbind(1, rep(0.5, length(L1pred)))##matrix(1, nrow = length(L1pred))
Z1pred <- Z2pred <- bbase(x = L1pred, xl = x1, xr = x2, ndx = 20, deg = 3)

## random effects
## hauls are the same for both levels
G1 <- G2 <- model.matrix( ~ -1 + factor(Haul), data = dat)

## fit in TMB
compile("nested_gam.cpp")
dyn.load(dynlib("nested_gam"))    
 
obj <- MakeADFun(
    data = list(Y1 = Y1,
                X1 = X1,
                Z1 = Z1,
                O1 = O1,
                G1 = G1,
                X1pred = X1pred,
                Z1pred = Z1pred,
                Y2 = Y2,
                X2 = X2,
                Z2 = Z2,
                O2 = O2,
                G2 = G2,
                X2pred = X2pred,
                Z2pred = Z2pred
                ),
    parameters = list(beta1 = matrix(rep(0, 2)),
                      lnsigmau1 = log(0.1),
                      u1 = matrix(0, nrow = ncol(Z1)),
                      lnsigmab1 = log(0.1),
                      b1 = matrix(0, nrow = ncol(G1)),
                      beta2 = matrix(rep(0, 2)),
                      lnsigmau2 = log(0.1),
                      u2 = matrix(0, nrow = ncol(Z2)),
                      lnsigmab2 = log(0.1),
                      b2 = matrix(0, nrow = ncol(G2))                      
                      ),
    random = c("u1", "u2", "b1", "b2"),
    DLL = "nested_gam"
)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = list(eval.max = 1e3, iter.max = 1e3))
rep <- sdreport(obj)
srep <- summary(rep)

## random effects plot
hauls <- as.numeric(gsub("factor\\(Haul\\)", "", colnames(G1)))
b1.hat <- srep[rownames(srep) == "b1", "Estimate"]
b2.hat <- srep[rownames(srep) == "b2", "Estimate"]
re.df <- rbind(data.frame(spp = spp, Level = "Test/Control", Haul = hauls, re = b1.hat),
               data.frame(spp = spp, Level = "Lower/Upper", Haul = hauls, re = b2.hat))
assign(paste0(spp, "_re.df"), re.df)

## CONTROL VS TEST
etaC <- as.data.frame(srep[rownames(srep) == "eta1pred", ])
etaC[, 1] <- - etaC[, 1]
rownames(etaC) <- NULL
names(etaC) <- c("eta", "eta.se")
etaC$Length <- L1pred
etaC$variable <- "Control"
etaT <- as.data.frame(srep[rownames(srep) == "eta1pred", ])
rownames(etaT) <- NULL
names(etaT) <- c("eta", "eta.se")
etaT$Length <- L1pred
etaT$variable <- "Test"
##
CTpredict <- rbind(etaC, etaT)
CTpredict$Proportion <- with(CTpredict, plogis(eta))
CTpredict$Lower <- with(CTpredict, plogis(eta + qnorm(0.025) * eta.se))
CTpredict$Upper <- with(CTpredict, plogis(eta + qnorm(0.975) * eta.se))
## work out marginals later
## for now working with conditionals
## 
dat$Control_Raised <- with(dat, Control_Count / Control_SUBSRATIO)
dat$Lower_Raised <- with(dat, Lower.codend_Count / Lower.codend_SUBSRATIO)
dat$Upper_Raised <- with(dat, Upper.codend_Count / Upper.codend_SUBSRATIO)
dat$Test_Raised <- with(dat, Lower_Raised + Upper_Raised)
dat$N <- with(dat, Control_Raised + Test_Raised)
dat$propC <- with(dat, Control_Raised / N)
dat$propT <- with(dat, Test_Raised / N)
library(reshape)
dat.melt <- melt(dat[, c("Haul", "Length", "propC", "propT")], id.var = c("Haul", "Length"))
names(dat.melt)[names(dat.melt) == "value"] <- "Proportion"
dat.melt$variable <- as.character(dat.melt$variable)
dat.melt$variable[dat.melt$variable == "propC"] <- "Control"
dat.melt$variable[dat.melt$variable == "propT"] <- "Test"
dat.melt$variable <- factor(dat.melt$variable, levels = c("Control", "Test"))
dat.melt$fHaul <- factor(dat.melt$Haul)
dat.melt <- merge(dat.melt, dat[, c("Haul", "Length", "N")])
## totals
totals <- aggregate(.~Length, FUN = sum,  data = dat[, c("Length", "Control_Raised", "Test_Raised")])
totals$N <- with(totals, Test_Raised + Control_Raised)
## proportions
totals$Control <- with(totals, Control_Raised / N)
totals$Test <- with(totals, Test_Raised / N)
totals.long <- melt(totals[, c("Length", "Test", "Control")], id.vars = "Length")
names(totals.long)[names(totals.long) == "value"] <- "Proportion"
## bring in N
totals.long <- merge(totals.long, totals[, c("Length", "N")])
## quick gam fit
dat$o <- with(dat, log(Test_SUBSRATIO / Control_SUBSRATIO))
gfit <- gam(N1 ~ s(Length, bs = "cr") + offset(o), data = dat, family = binomial, method = "ML")
pred.df <- data.frame(Length = L1pred, o = 0)
pred.df$Test <- predict(gfit, newdata = pred.df, type = "response")
pred.df$Control <- 1 - pred.df$Test
pred.long <- melt(pred.df[, names(pred.df) != "o"], id.vars = "Length")
names(pred.long)[names(pred.long) == "value"] <- "Proportion"

xlab.val <- ifelse(spp == "Nephrops", "Carapace length (mm)", "Length (cm)")

p <- ggplot(dat.melt, aes(x = Length, y = Proportion)) +
    geom_point(aes(colour = fHaul, size = log(N)), alpha = 0.6) +
    facet_wrap(~ variable, ncol = 2) +
    theme(legend.position = "none") +
    geom_ribbon(data = CTpredict, aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
    geom_line(data = CTpredict, colour = "navy") +
    geom_hline(yintercept = 0.5, colour = "grey") +
    ##geom_line(data = pred.long, colour = "navy", linetype = 2) +
    geom_point(data = totals.long, aes(size = log(N)), pch = 21, fill = "gold") +
    ##geom_point(data = totals.long, pch = 21, fill = "gold") +
    scale_size(range = c(0, 4)) +
    xlab(xlab.val)## +
    ##ggtitle(spp)

##
pdf(paste0("../tex/figures/", spp, "_control-test_nested_v2.pdf"), height = 6, width = 12)
print(p)
dev.off()

## LOWER VS UPPER
etaL <- as.data.frame(srep[rownames(srep) == "eta2pred", ])
rownames(etaL) <- NULL
names(etaL) <- c("eta", "eta.se")
etaL$Length <- L1pred
etaL$variable <- "Lower"
etaU <- as.data.frame(srep[rownames(srep) == "eta2pred", ])
etaU[,1] <- -etaU[,1]
rownames(etaU) <- NULL
names(etaU) <- c("eta", "eta.se")
etaU$Length <- L1pred
etaU$variable <- "Upper"
LUpredict <- rbind(etaL, etaU)
LUpredict$Proportion <- with(LUpredict, plogis(eta))
LUpredict$Lower <- with(LUpredict, plogis(eta + qnorm(0.025) * eta.se))
LUpredict$Upper <- with(LUpredict, plogis(eta + qnorm(0.975) * eta.se))
LUpredict$variable <- factor(LUpredict$variable, levels = c("Upper", "Lower"))
## 
library(reshape)
dat$propL <- with(dat, Lower_Raised / (Lower_Raised + Upper_Raised))
dat$propU <- with(dat, Upper_Raised / (Lower_Raised + Upper_Raised))
dat$N <- with(dat, Lower_Raised + Upper_Raised)
dat.melt <- melt(dat[, c("Haul", "Length", "propL", "propU")], id.var = c("Haul", "Length"))
names(dat.melt)[names(dat.melt) == "value"] <- "Proportion"
dat.melt$variable <- as.character(dat.melt$variable)
dat.melt$variable[dat.melt$variable == "propL"] <- "Lower"
dat.melt$variable[dat.melt$variable == "propU"] <- "Upper"
dat.melt$variable <- factor(dat.melt$variable, levels = c("Lower", "Upper"))
dat.melt$fHaul <- factor(dat.melt$Haul)
dat.melt <- merge(dat.melt, dat[, c("Haul", "Length", "N")])
dat.melt$variable <- relevel(dat.melt$variable, ref = "Upper")

## totals
totals <- aggregate(.~Length, FUN = sum,  data = dat[, c("Length", "Lower_Raised", "Upper_Raised")])
totals$N <- with(totals, Upper_Raised + Lower_Raised)
## proportions
totals$Lower <- with(totals, Lower_Raised / N)
totals$Upper <- with(totals, Upper_Raised / N)
totals.long <- melt(totals[, c("Length", "Upper", "Lower")], id.vars = "Length")
names(totals.long)[names(totals.long) == "value"] <- "Proportion"
## bring in N
totals.long <- merge(totals.long, totals[, c("Length", "N")])

## quick gam fit
dat$o <- with(dat, log(Lower.codend_SUBSRATIO / Upper.codend_SUBSRATIO))
gfit <- gam(N2 ~ s(Length, bs = "cr") + offset(o), data = dat, family = binomial, method = "ML")
pred.df <- data.frame(Length = L1pred, o = 0)
pred.df$Lower <- predict(gfit, newdata = pred.df, type = "response")
pred.df$Upper <- 1 - pred.df$Lower
pred.long <- melt(pred.df[, names(pred.df) != "o"], id.vars = "Length")
names(pred.long)[names(pred.long) == "value"] <- "Proportion"
pred.long$variable <- relevel(pred.long$variable, ref = "Upper")

p <- ggplot(dat.melt, aes(x = Length, y = Proportion)) +
    geom_point(aes(colour = fHaul, size = log(N)), alpha = 0.6) +
    facet_wrap(~ variable, ncol = 1) +
    theme(legend.position = "none") +
    geom_ribbon(data = LUpredict, aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
    geom_line(data = LUpredict, colour = "navy") +
    ##geom_line(data = pred.long, colour = "navy", linetype = 2) +
    geom_point(data = totals.long, aes(size = log(N)), pch = 21, fill = "gold") +
    geom_hline(yintercept = 0.5, , colour = "grey") +
    scale_size(range = c(0, 4)) +
    xlab(xlab.val)
    ##ggtitle(spp)

pdf(paste0("../tex/figures/", spp, "_lower-upper_nested_v2.pdf"), height = 6, width = 6)
print(p)
dev.off()
##}  

re.df <- rbind(Haddock_re.df, Nephrops_re.df, Whiting_re.df)

pdf("../tex/figures/dual_codend_nested_random_effects.pdf", height = 9, width = 6)
ggplot(re.df, aes(x = Haul, y = re)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(shape = Level, fill = Level), size = 2) +
    facet_wrap(~spp, ncol = 1, scales = "free_y") +
    scale_shape_manual(values = c(21, 17)) +
    scale_x_continuous(breaks = hauls) +
    scale_fill_manual(values = c("white", "black")) +
    theme(legend.position = "bottom",
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    ylab("Random effect")
dev.off()
