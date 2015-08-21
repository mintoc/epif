library("MASS")
library("gdata")
library("reshape") 
library("ggplot2")
library("MGLM")

##Load in data
nw.dat<- read.xlsx("C:\\Users\\bburke\\Desktop\\GMIT\\Gear Trials\\Data\\2015 BIM Nephrops quad rig trials\\Arkh Angell 70 V 100mm codends the Smalls 2015\\Nephrops.xlsx", sheetName="Nephrops weight", header=TRUE, perl="C:\\Strawberry\\perl\\bin\\perl.exe")
nl.dat<- read.xlsx("C:\\Users\\bburke\\Desktop\\GMIT\\Gear Trials\\Data\\2015 BIM Nephrops quad rig trials\\Arkh Angell 70 V 100mm codends the Smalls 2015\\Nephrops.xlsx", sheetName="Nephrops length", header=TRUE, perl="C:\\Strawberry\\perl\\bin\\perl.exe")

##Some data cleaning
nl.dat<-subset(nl.dat, select=c("Haul", "Compartment", "Gear.no.", "Species", "Carapace.length", "Count", "Raising.factor", "Raised.count"))
nw.dat<-subset(nw.dat, select=-perl)
nl.dat$Compartment <- gsub(" ", "", nl.dat$Compartment)
nw.dat$Compartment <- gsub(" ", "", nw.dat$Compartment)
nl.dat$Haul<-paste("H", nl.dat$Haul, sep="")
nw.dat$Haul<-paste("H", nw.dat$Haul, sep="")

##Merge the two sheets together
nlnw<- merge(nl.dat, nw.dat, by=c("Haul", "Compartment"))
nlnw<-subset(nlnw, select=-Species.y)
nlnw<-subset(nlnw, select=c("Haul", "Compartment", "Carapace.length", "Count", "Raising.factor.x", "Sample.wt"))

##Remove the 80mm mesh size hauls
##Hauls 11 and 12 can be removed as there was no 100mm tow associated with them
nlnw<-subset(nlnw, select=Compartment!="80mmPort")
nlnw<-subset(nlnw, select=Compartment!="80mmStarboard")
nlnw<-subset(nlnw, !Haul %in% c("H11", "H12"))

##Next prepare the Subsratio
nlnw$Raising.factor.x[nlnw$Haul == "H4" & 
                        nlnw$Compartment == "100mmPort"] <- 6.8
nlnw$Raising.factor.x[nlnw$Haul == "H9" & 
                        nlnw$Compartment == "100mmStarboard"] <- 15.4

nlnw$SUBSRATIO<-1/nlnw$Raising.factor.x

##Reshaping data frame
vars2keep<- c("Haul", "Compartment", "Carapace.length", "Count", "SUBSRATIO")
AA.melt<-melt(nlnw[,vars2keep], id=c("Compartment", "Carapace.length", "Haul") )
AA.cast<-cast(AA.melt, Haul + Carapace.length ~ Compartment+variable)

##Rename Columns
names(AA.cast)[names(AA.cast)=="100mmPort_Count"]<-"c100mm_Port"
names(AA.cast)[names(AA.cast)=="100mmStarboard_Count"]<-"c100mm_Starboard"
names(AA.cast)[names(AA.cast)=="70mmPort_Count"]<-"c70mm_Port"
names(AA.cast)[names(AA.cast)=="70mmStarboard_Count"]<-"c70mm_Starboard"

names(AA.cast)[names(AA.cast)=="100mmPort_SUBSRATIO"]<-"s100mm_Port"
names(AA.cast)[names(AA.cast)=="100mmStarboard_SUBSRATIO"]<-"s100mm_Starboard"
names(AA.cast)[names(AA.cast)=="70mmPort_SUBSRATIO"]<-"s70mm_Port"
names(AA.cast)[names(AA.cast)=="70mmStarboard_SUBSRATIO"]<-"s70mm_Starboard"

##Create Variables for Carapace Length squared and cubed
AA.cast$pCarapace.length<-AA.cast$Carapace.length/max(AA.cast$Carapace.length)
AA.cast$pCarapace.length2<-AA.cast$pCarapace.length^2
AA.cast$pCarapace.length3<-AA.cast$pCarapace.length^3

##Create loop to address the missing Subsratio
for(i in 1:dim(AA.cast)[1]){
  haul.dat <- subset(AA.cast, Haul == AA.cast$Haul[i])
  
  if(is.na(AA.cast$s100mm_Port[i])){
    AA.cast$s100mm_Port[i] <- unique(na.omit(haul.dat$s100mm_Port))
  }
  
  if(is.na(AA.cast$s100mm_Starboard[i])){
    AA.cast$s100mm_Starboard[i] <- unique(na.omit(haul.dat$s100mm_Starboard))
  }
  
  if(is.na(AA.cast$s70mm_Port[i])){
    AA.cast$s70mm_Port[i] <- unique(na.omit(haul.dat$s70mm_Port))
  }
  
  if(is.na(AA.cast$s70mm_Starboard[i])){
    AA.cast$s70mm_Starboard[i] <- unique(na.omit(haul.dat$s70mm_Starboard))
  }
}  

AA.cast[is.na(AA.cast)]<-0

##incorporate subsratio
##create a subsratio matrix - ref notebook on layout 
subs.var<- c("s100mm_Port", "s100mm_Starboard", "s70mm_Port", "s70mm_Starboard")

sub.mat<-as.matrix(AA.cast[,subs.var])

colnames(sub.mat)<-c("s100mm_Port", "s100mm_Starboard", "s70mm_Port", "s70mm_Starboard")

ggt <- cbind(AA.cast$s100mm_Port, AA.cast$s100mm_Starboard, AA.cast$s70mm_Port, AA.cast$s70mm_Starboard)

##Set up offset
offset.mat <- log(apply(ggt, MARGIN = 2, FUN = function(zz){zz/ggt[,1]}))

##Try to use AA.cast as the dataset for the Dirichlet multinomial
reg<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
          pCarapace.length, data=AA.cast, dist="DM")

reg2<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                pCarapace.length+pCarapace.length2, data=AA.cast, dist="DM")

reg2.mnom<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                 pCarapace.length+pCarapace.length2, data=AA.cast, dist="MN")

reg3<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                pCarapace.length+pCarapace.length2+pCarapace.length3, data=AA.cast, dist="DM")

reg4<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                 pCarapace.length+pCarapace.length2+pCarapace.length3+offset(offset.mat), data=AA.cast, dist="DM")

reg4<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                 1 + offset(offset.mat), data=AA.cast, dist="DM")

reg4.no.offset<- MGLMreg(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                 1 , data=AA.cast, dist="DM")

reg4<- multinom(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                 1 + offset(offset.mat), data=AA.cast)

reg4.no.offset<- multinom(cbind(c100mm_Port, c100mm_Starboard, c70mm_Port, c70mm_Starboard)~
                           1 , data=AA.cast)


mnom.fit <- multinom(y ~ x + offset(offset.mat))

AIC.models<- c(reg$AIC,reg2$AIC,reg3$AIC, reg4$AIC)

##Plot Dirichlet Model for reg, reg2 and reg3
pcl<-10:59/max(nlnw$Carapace.length)

predX<-cbind(1, pcl)
predY<-cbind(1, pcl, pcl^2)
predZ<-cbind(1, pcl, pcl^2, pcl^3)

pred.prop2 <- predict(reg, newdata = predX)
pred.prop3 <- predict(reg2, newdata = predY)
pred.prop4 <- predict(reg3, newdata = predZ)

nn <- dim(pred.prop2)[1]
tt <- dim(pred.prop3)[1]
ll <- dim(pred.prop4)[1]

pred.DM.df2<-data.frame(
  Mesh.Size=factor(rep(
    c("c100mmPort", "c100mmStarboard", "c70mmPort", "c70mmStarboard"), 
    each=nn), levels=c("c100mmPort", "c100mmStarboard", "c70mmPort", "c70mmStarboard")),
  Carapace.length=rep(10:59, times=4),
  proportion=c(pred.prop2)
)

plot + geom_line(data = pred.DM.df2)

pred.DM.df3<-data.frame(
  Mesh.Size=factor(rep(
    c("c100mmPort", "c100mmStarboard", "c70mmPort", "c70mmStarboard"), 
    each=tt), levels=c("c100mmPort", "c100mmStarboard", "c70mmPort", "c70mmStarboard")),
  Carapace.length=rep(10:59, times=4),
  proportion=c(pred.prop3)
)

plot + geom_line(data = pred.DM.df3)

pred.DM.df4<-data.frame(
  Mesh.Size=factor(rep(
    c("c100mmPort", "c100mmStarboard", "c70mmPort", "c70mmStarboard"), 
    each=ll), levels=c("c100mmPort", "c100mmStarboard", "c70mmPort", "c70mmStarboard")),
  Carapace.length=rep(10:59, times=4),
  proportion=c(pred.prop4)
)

plot + geom_line(data = pred.DM.df4)




