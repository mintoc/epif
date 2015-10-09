

bulk.vec <- unlist(bulk.cast[, -1])
##cuts <- c(0, quantile(bulk.vec, c(1/3,2/3)), 2 * max(bulk.vec))
cuts <- c(0, quantile(bulk.vec, c(0.5)), 2 * max(bulk.vec))

unique(as.data.frame(apply(bulk.cast[,-1], 2, cut, breaks = cuts)))
