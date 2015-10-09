
## function to plot residuals
plot.multinom.resid <- function(covar, ylim = NULL){
  par(mfrow = c(2,2), mar = c(2,2,1,1), oma = c(2,2,1,1))
  x.vec <- our.lass.neph.cast[, covar] 
  for(i in 1:4){
    if(class(x.vec) %in% c("numeric", "integer")){
      if(is.null(ylim)){
        plot(x.vec, mnom.resid[,i], pch = 19, col = "grey", cex = log10(neph.count.mat[,i] + 1.2))
      }else{
        plot(x.vec, mnom.resid[,i], pch = 19, col = "grey", cex = log10(neph.count.mat[,i] + 1.2), ylim = ylim)
      }
      ##lo <- loess(mnom.resid[,i] ~ x.vec, weights = neph.count.mat[,i], span = 0.75)
      lo <- loess(mnom.resid[,i] ~ x.vec, span = 0.75)
      curve(predict(lo, newdata = data.frame(x.vec = x)), add = TRUE, col = "red", lwd = 2)
    }
    if(class(x.vec) == "factor"){
      if(is.null(ylim)){
        boxplot(mnom.resid[,i] ~ x.vec, notch = TRUE, col = "grey")
      }else{
        boxplot(mnom.resid[,i] ~ x.vec, notch = TRUE, col = "grey", ylim = ylim)
      }
    }
    abline(h = 0)
    title(colnames(neph.count.mat)[i])
  }
}
