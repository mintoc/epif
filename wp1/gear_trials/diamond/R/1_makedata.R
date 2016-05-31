##-------------------
## Quad-rig data
## CM, BB: Tue Sep 22 2015
## Note: bigger models run on the cluster
##-------------------
library(gdata)

##----------
## OUR LASS 
##----------
our.lass.neph.dat <- read.xls("../data/2015 BIM Nephrops quad rig trials/Our Lass 2 70_80_90_100mm codends Irish Sea July 2015/Nephrops Raised Counts Our Lass 2 Irish Sea July 2015.xlsx", 
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

bulk.cast.diff <- bulk.cast[,-1] - bulk.cast[,"m70mm_Total.catch"]
names(bulk.cast.diff) <- paste("d", names(bulk.cast.diff), sep = "")
bulk.cast.diff$fHAUL <- bulk.cast$fHAUL

## merge counts and subsampling ratio back together 
our.lass.neph.cast0 <- merge(our.lass.neph.cast, subs.cast, by = "fHAUL", all.x = TRUE)

our.lass.neph.cast1 <- merge(our.lass.neph.cast0, bulk.cast, by = "fHAUL", all.x = TRUE)

our.lass.neph.cast <- merge(our.lass.neph.cast1, netpos.cast, by = "fHAUL", all.x = TRUE)

our.lass.neph.cast <- our.lass.neph.cast[order(our.lass.neph.cast$fHAUL, our.lass.neph.cast$Carapace.length),]

## show first few lines
head(our.lass.neph.cast, 2)

## Create "net configuration" variable
nc.vars <- c("m70mm_Net.position", "m80mm_Net.position", "m90mm_Net.position", "m100mm_Net.position")
our.lass.neph.cast$netconfig <- factor(paste("NC", apply(our.lass.neph.cast[, nc.vars], 1, paste, collapse = ""), sep =""))

###Position matrices

net.names <- c("m70mm_Net.position", "m80mm_Net.position", "m90mm_Net.position", "m100mm_Net.position")

PO <- (our.lass.neph.cast[,net.names] == 1) * 1
PI <- (our.lass.neph.cast[,net.names] == 2) * 1
SI <- (our.lass.neph.cast[,net.names] == 3) * 1
SO <- (our.lass.neph.cast[,net.names] == 4) * 1

colnames(PO) <- colnames(PI) <- colnames(SI) <- colnames(SO) <- NULL
rownames(PO) <- rownames(PI) <- rownames(SI) <- rownames(SO) <- NULL

## Extract the matrix of counts
count.vars <- c("m70mm_Count", "m80mm_Count", "m90mm_Count", "m100mm_Count")

neph.count.mat <- as.matrix(our.lass.neph.cast[, count.vars])

## Extract the matrix of subsampling ratios
subsratio.vars <- c("m70mm_SUBSRATIO", "m80mm_SUBSRATIO", "m90mm_SUBSRATIO", "m100mm_SUBSRATIO")

subsratio.mat <- as.matrix(our.lass.neph.cast[, subsratio.vars])

## Create the offset (NEED TO CHECK THIS)
offset.mat <- log(apply(subsratio.mat, 2, FUN = 
                          function(zz){zz/subsratio.mat[,1]}))

## SAVE OBJECTS
save(list = c("our.lass.neph.cast", "offset.mat", "neph.count.mat", "count.vars", "subsratio.vars", "nc.vars", "PO", "PI", "SI", "SO"), file = "our_lass_data_objects.RData")
