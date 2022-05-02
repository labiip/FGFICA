GetPatterns <- function(ICA_result, confidence.var) {
  library(PDFEstimator)
  
  analysis.result <- as.matrix(ICA_result)
  up.index <- data.frame()
  bottom.index <- data.frame()
  up.var <- 1
  bottom.var <- 1
  for (variable in 1:ncol(analysis.result)) {
    cat("Col ID:", variable, "\t")
    temp <- analysis.result[, variable]
    
    dist = estimatePDF(temp)
    UpLimIndex <- min(which(dist$cdf >= (1 - (1 - confidence.var)/2)))
    UpLim <- dist$x[UpLimIndex]
    cat("Up Limit:", UpLim, "\t")
    BottomLimIndex <- max(which(dist$cdf <= ((1 - confidence.var)/2)))
    BottomLim <- dist$x[BottomLimIndex]
    cat("Battom Limit:", BottomLim, "\t")
    
    temp.up <- which(temp >= UpLim)
    if (length(temp.up) != 0) {
      cat("Up length:", length(temp.up), "\t")
      up.index[1:length(temp.up), up.var] <- temp.up
      up.var <- up.var + 1
    }
    
    temp.bottom <- which(temp <= BottomLim)
    if (length(temp.bottom) != 0) {
      cat("bottom length:", length(temp.bottom), "\n")
      bottom.index[1:length(temp.bottom), bottom.var] <- temp.bottom
      bottom.var <- bottom.var + 1
    }
  }
  up.index <- as.matrix(up.index)
  bottom.index <- as.matrix(bottom.index)
  return(list(up.index = up.index, bottom.index = bottom.index))
}
