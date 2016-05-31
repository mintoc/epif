###JOHN HINDE'S MODELS FOR THE MINERS DATA

##load("C:/Users/bburke/Desktop/GMIT/R/Statistics in R - J Hinde/SMIR/data/R_data/miners.rda")
load("/media/sf_docs/analyses/epif/wp1/gear_trials/diamond/R/miners/SMIR/data/miners.rda")
miners <- transform(miners, t = n + m + s)
miners <- transform(miners, np=n/t, mp=m/t, sp=s/t)


###################################################
### chunk number 20:  eval=FALSE
###################################################
## data(miners)
## miners <- transform(miners, t = n + m + s)
## miners <- transform(miners, np=n/t, mp=m/t, sp=s/t)


###################################################
### chunk number 21: miners.plot eval=FALSE
###################################################
## miners.plot <- xyplot(np+mp+sp~years,data=miners,col=1,
## ylab="proportion",type="p",pch=c("n","m","s"))
## print(update(miners.plot,type="b",cex=1.5,lty=1:3))
## #print(xyplot(np+mp+sp~years,data=miners,col=1,cex=2,
## #ylab="proportion",type="b",pch=c("n","m","s"),lty=c(1,2,3)))


###################################################
### chunk number 22: 
###################################################
miners.plot <- xyplot(np+mp+sp~years,data=miners,col=1,
                      ylab="proportion",type="p",pch=c("n","m","s"))

print(update(miners.plot,type="b",cex=1.5,lty=1:3))
#print(xyplot(np+mp+sp~years,data=miners,col=1,cex=2,
#ylab="proportion",type="b",pch=c("n","m","s"),lty=c(1,2,3)))


###################################################
### chunk number 23: miners.glm
###################################################
miners <- transform(miners,period=factor(years))
miners.glm <- glm(cbind(n, (t - n)) ~ period, family = binomial, 
                  data = miners)
summary(miners.glm)


###################################################
### chunk number 24: miners.glm1
###################################################
miners.glm1 <- update(miners.glm,. ~ years)
#summary(miners.glm1)
print(anova(miners.glm1),digits=4)
round(resid(miners.glm1,type='pearson'),3)


###################################################
### chunk number 25: miners.glm1a
###################################################
miners.glm1a <- update(miners.glm1, family = binomial(link = "cloglog"))
#summary(miners.glm1a)
print(anova(miners.glm1a),digits=4)
round(resid(miners.glm1a,type='pearson'),3)


###################################################
### chunk number 26: miners.glm2
###################################################
miners <- transform(miners, ly=log(years))
miners.glm2 <- update(miners.glm1a,.~ly)
summary(miners.glm2)
round(resid(miners.glm2,type="pearson"),3)


###################################################
### chunk number 27: 
###################################################
miners$cnp <- fitted(miners.glm2)


###################################################
### chunk number 28: miners.glm3
###################################################
miners.glm3 <- update(miners.glm2,.~ly,family=binomial)
summary(miners.glm3)
round(resid(miners.glm3,type='pearson'),3)


###################################################
### chunk number 29: miners.plot2 eval=FALSE
###################################################
## miners$lnp <- fitted(miners.glm3)
## print(xyplot(np+cnp+lnp~years,type=c("p","l","l"),data=miners,ylab="proportion",
## panel=panel.superpose.2,lty=c(1,2)))


###################################################
### chunk number 30: miners.plot3 eval=FALSE
###################################################
## print(xyplot(np~years,data=miners,ylab='proportion',ylim=c(0,1.2),col="black",
## panel=function(x,y){
## panel.xyplot(x,y)
## panel.curve(1-exp(-exp(4.582-1.312*log(x))))
## panel.curve((eta <- exp(9.609-2.576*log(x)))/(1+eta),lty=2)
## }))


###################################################
### chunk number 31: 
###################################################
print(xyplot(np~years,data=miners,ylab='proportion',ylim=c(0,1.2),col="black",
             panel=function(x,y){
               panel.xyplot(x,y)
               panel.curve(1-exp(-exp(4.582-1.312*log(x))))
               panel.curve((eta <- exp(9.609-2.576*log(x)))/(1+eta),lty=2)
             }))


###################################################
### chunk number 32: miners.empirical.logits
###################################################
miners <- transform(miners,logit2 = log(m/n), logit3 = log(s/n))
print(miners[c("years","logit2","logit3")])


###################################################
### chunk number 33: miners.logit.plot eval=FALSE
###################################################
## print(xyplot(logit2+logit3~years,type="b",pch=c("m","s"),
## data=miners,cex=1.5, ylab="logit",lty=c(1,2),lwd=1.8))


###################################################
### chunk number 34: 
###################################################
miners$lnp <- fitted(miners.glm3)
print(xyplot(np+cnp+lnp~years,type=c("p","l","l"),data=miners,ylab="proportion",
             panel=panel.superpose.2,lty=c(1,2)))


###################################################
### chunk number 35: 
###################################################
print(xyplot(logit2+logit3~years,type="b",pch=c("m","s"),
             data=miners,cex=1.5, ylab="logit",lty=c(1,2),lwd=1.8))


###################################################
### chunk number 36: miners.log.plot eval=FALSE
###################################################
## print(xyplot(logit2+logit3~years,type="b",pch=c("m","s"),col="black",lty=c(2,3),lwd=1.5,
## data=miners,cex=1.5,ylab="logit",xlab="years (log scale)",scale=list(x=list(log=TRUE,
## at=c(10,20,30,40,50),labels=c("10","20","30","40","50")))))


###################################################
### chunk number 37: 
###################################################
print(xyplot(logit2+logit3~years,type="b",pch=c("m","s"),col="black",lty=c(2,3),lwd=1.5,
             data=miners,cex=1.5,ylab="logit",xlab="years (log scale)",scale=list(x=list(log=TRUE,
                                                                                         at=c(10,20,30,40,50),labels=c("10","20","30","40","50")))))
###################################################
### chunk number 42: miners.multinomial
###################################################
miners.long <- reshape(miners,drop=names(miners)[6:11],
                       varying=list(names(miners)[2:4]),
                       timevar="resp", times=c("n","m","s"),
                       v.names="count", direction="long")
miners.long <- transform(miners.long, resp=factor(resp))
miners.long <- transform(miners.long, resp = relevel(resp,ref="n"),
                         lyr = log(years),period = factor(years))
miners.glm5 <- glm(count~period+resp,data=miners.long,family=poisson)

#miners.gnm1 <- gnm(count~resp,eliminate=period,data=miners.long,family=poisson,ofInterest="resp")
coef(summary(miners.glm5))[9:10,,drop=FALSE]
#summary(miners.gnm1)
miners.glm6 <- update(miners.glm5,.~.+resp:years)
#miners.gnm2 <- update(miners.gnm1,.~.+resp:years)
anova(miners.glm5,miners.glm6)
#anova(miners.gnm1,miners.gnm2)
round(summary(miners.glm6)$coef[9:12,],4)
#round(coef(summary(miners.glm6))[9:12,,drop=FALSE],4)
#summary(miners.gnm2)
#getContrasts(miners.gnm2,ofInterest(miners.gnm2)[5:3])


###################################################
### chunk number 43: miners.glm7
###################################################
miners.glm7 <- update(miners.glm6,.~.+years+resp:years)
#miners.gnm3 <- update(miners.gnm2,.~.+years+resp:years)
print(anova(miners.glm7,miners.glm6),digits=5)
#anova(miners.6nm2,miners.gnm3)
round(coef(summary(miners.glm7))[9:12,],4)
#summary(miners.gnm3)
miners.long$fval7 <- fitted(miners.glm7)
miners.long$resid7 <- resid(miners.glm7)
options(digits=3)
cbind(miners.long[1:8,c('count','fval7','resid7')],
      miners.long[9:16,c('count','fval7','resid7')],
      miners.long[17:24,c('count','fval7','resid7')],
      row.names=NULL)


###################################################
### chunk number 44: miners.glm8
###################################################
miners.glm8 <- update(miners.glm6,.~period+resp+lyr+resp:lyr)
#miners.gnm4 <- update(miners.gnm3,.~resp*lyr)
print(anova(miners.glm6,miners.glm8),digits=4)
#anova(miners.gnm2,miners.gnm4)
print(coef(summary(miners.glm8))[9:12,,drop=FALSE],digits=5)
#summary(miners.gnm4)
miners.long$fval8 <- fitted(miners.glm8)
miners.long$resid8 <- resid(miners.glm8)
options(digits=3)
cbind(miners.long[1:8,c('count','fval8','resid8')],
      miners.long[9:16,c('count','fval8','resid8')],
      miners.long[17:24,c('count','fval8','resid8')],row.names=NULL)


###################################################
### chunk number 45: miners.plot5 eval=FALSE
###################################################
## miners.long$fv8 <- fitted(miners.glm8)
## #miners.long$fv <- fitted(miners.gnm4)
## miners.long <- transform(miners.long,
##                          f2p8=fv8/t,
##                          op = count/t)
## print(miners.plot <- xyplot(op~years,data=miners.long,groups=resp,
## ylab='proportion',xlab='years',pch=c("n","m","s"),cex=1.5))
## print(update(miners.plot,
## panel=function(x,y,...){
## panel.xyplot(x,y,...)
## panel.xyplot(x,miners.long$f2p8,type="l",...)}))


###################################################
### chunk number 46: 
###################################################
miners.long$fv8 <- fitted(miners.glm8)
#miners.long$fv <- fitted(miners.gnm4)
miners.long <- transform(miners.long,
                         f2p8=fv8/t,
                         op = count/t)
print(miners.plot <- xyplot(op~years,data=miners.long,groups=resp,
                            ylab='proportion',xlab='years',pch=c("n","m","s"),cex=1.5))
print(update(miners.plot,
             panel=function(x,y,...){
               panel.xyplot(x,y,...)
               panel.xyplot(x,miners.long$f2p8,type="l",...)}))


###################################################
### chunk number 47: 
###################################################
miners.long$respm <- relevel(miners.long$resp,ref="m")
coef(summary(update(miners.glm8,.~period+respm*lyr)))[12,,drop=FALSE]
#getContrasts(miners.gnm4,ofInterest(miners.gnm4)[3:4])


###################################################
### chunk number 48: miners.glm9
###################################################
miners.long <- transform(miners.long,r23=ifelse(resp!="n",1,0))
miners.long <- transform(miners.long,r23lyr=r23*lyr)
miners.glm9 <- update(miners.glm8,.~factor(years)+resp+r23lyr)
#miners.gnm5 <- update(miners.gnm4,.~resp+r23:lyr)
anova(miners.glm9,miners.glm8)
#anova(miners.gnm5,miners.gnm4)
print(coef(summary(miners.glm9))[9:11,,drop=FALSE],5)
#coef(summary(miners.gnm5))[-(1:8),]


###################################################
### chunk number 49: miners.plot6 eval=FALSE
###################################################
## miners.long$fval9 <- fitted(miners.glm9)
## miners.long <- transform(miners.long,
##                          f2p9=fval9/t,
##                          op = count/t)
## print(update(miners.plot,
## panel=function(x,y,...){
## panel.xyplot(x,y,...)
## panel.xyplot(x,miners.long$f2p9,type="l",...)}))


###################################################
### chunk number 50: 
###################################################
miners.long$fval9 <- fitted(miners.glm9)
miners.long <- transform(miners.long,
                         f2p9=fval9/t,
                         op = count/t)
print(update(miners.plot,
             panel=function(x,y,...){
               panel.xyplot(x,y,...)
               panel.xyplot(x,miners.long$f2p9,type="l",...)}))
