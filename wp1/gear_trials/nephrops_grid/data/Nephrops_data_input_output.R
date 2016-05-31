##-------------------------------------------
## read in and output the Nephrops data in required format
## Thu Mar 24 2016
## Notes: please check a few examples to make sure it's correct
##-------------------------------------------

neph.dat <- read.csv("Neph length 24_03_2016.csv", header = FALSE)

## first row is haul number
## second row is compartment number

data.index <- 3:nrow(neph.dat)

## example: for the first column
data.frame(HAUL = neph.dat[1,1], COMPARTMENT = neph.dat[2,1], CARAPACE.LENGTH = neph.dat[data.index, 1])

## container data frame
neph.dat.format <- data.frame(NULL)

## loop through all
for(i in 1:ncol(neph.dat)){
    print(i)
    ## temporary data frame for the column i
    tmp <- data.frame(HAUL = neph.dat[1,i], COMPARTMENT = neph.dat[2,i], CARAPACE.LENGTH = neph.dat[data.index, i])
    ## remove zeros
    tmp <- subset(tmp, CARAPACE.LENGTH > 0)
    ## include in the overall data frame
    neph.dat.format <- rbind(neph.dat.format, tmp)
    ## clean up
    rm(tmp)
}

write.csv(neph.dat.format, file = "Neph length formatted 24_03_2016.csv", row.names = FALSE)

