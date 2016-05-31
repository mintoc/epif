library(gdata)
library(reshape2)

##our.lass.neph.dat <- read.xls("C:/Users/bburke/Desktop/GMIT/Gear Trials/Data/2015 BIM Nephrops quad rig trials/Our Lass 2 70_80_90_100mm codends Irish Sea July 2015/Nephrops Raised Counts Our Lass 2 Irish Sea July 2015.xlsx", 
##                              sheet = "All hauls", stringsAsFactors = FALSE)

our.lass.neph.dat <- read.xls("../data/2015 BIM Nephrops quad rig trials/Our Lass 2 70_80_90_100mm codends Irish Sea July 2015/Nephrops Raised Counts Our Lass 2 Irish Sea July 2015.xlsx", sheet = "All hauls", stringsAsFactors = FALSE)

our.lass.neph.dat <- subset(our.lass.neph.dat, Haul.No != 14)

net.pos.dat<- data.frame(haul = as.factor(our.lass.neph.dat$Haul.No),
                         mesh = our.lass.neph.dat$Mesh.Size, 
                         position = our.lass.neph.dat$Net.position)

data <- model.matrix(~mesh+position, net.pos.dat)
data <- cbind(net.pos.dat$haul, data)
colnames(data) <- c("haul", "intercept", "mesh70mm", "mesh80mm", "mesh90mm", "position")

po.mat<- subset(data, data[, "position"]==1)
rownames(po.mat) <- NULL
po.mat<-as.data.frame(po.mat)

pi.mat<- as.data.frame(subset(data, data[, "position"]==2))
rownames(pi.mat) <- NULL
pi.mat<-as.data.frame(pi.mat)

si.mat<- as.data.frame(subset(data, data[, "position"]==3))
rownames(si.mat) <- NULL
si.mat<-as.data.frame(si.mat)

so.mat<- as.data.frame(subset(data, data[, "position"]==4))
rownames(so.mat) <- NULL
so.mat<-as.data.frame(so.mat)

po.mat$mesh100mm<-1-rowSums(po.mat[,2:4])
pi.mat$mesh100mm<-1-rowSums(pi.mat[,2:4])
si.mat$mesh100mm<-1-rowSums(si.mat[,2:4])
so.mat$mesh100mm<-1-rowSums(so.mat[,2:4])

po.mat$intercept <- NULL
po.mat$position <- NULL

pi.mat$intercept <- NULL
pi.mat$position <- NULL

si.mat$intercept <- NULL
si.mat$position <- NULL

so.mat$intercept <- NULL
so.mat$position <- NULL

## same length

PO <- PI <- SI <- SO <- matrix(0, nrow = nrow(our.lass.neph.dat), ncol = 4,
                               dimnames = list(NULL, c("70mm", "80mm", "90mm", "100mm")))


## PORT OUTSIDE
for(i in 1:nrow(PO)){
  if(our.lass.neph.dat$Net.position[i] == 1){
    PO[i, our.lass.neph.dat$Mesh.Size[i]] <- 1
  }
}

## PORT INSIDE
for(i in 1:nrow(PI)){
  if(our.lass.neph.dat$Net.position[i] == 2){
    PI[i, our.lass.neph.dat$Mesh.Size[i]] <- 1
  }
}

## STARBOARD INSIDE
for(i in 1:nrow(SI)){
  if(our.lass.neph.dat$Net.position[i] == 3){
    SI[i, our.lass.neph.dat$Mesh.Size[i]] <- 1
  }
}

## STARBOARD OUTSIDE
for(i in 1:nrow(SO)){
  if(our.lass.neph.dat$Net.position[i] == 4){
    SO[i, our.lass.neph.dat$Mesh.Size[i]] <- 1
  }
}


## SANDBOX
## fast but unintuitive way of calculating the position matrices
po <- matrix(0, nrow = nrow(our.lass.neph.dat), ncol = 4,
             dimnames = list(NULL, c("70mm", "80mm", "90mm", "100mm")))

po.idx <- which(our.lass.neph.dat$Net.position == 1)

po.idx.mat <- cbind(po.idx,
                    match(our.lass.neph.dat[po.idx, "Mesh.Size"], colnames(po)))

po[po.idx.mat] <- 1
