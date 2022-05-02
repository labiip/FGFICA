FGFICA <- function(InputData, 
                   FusionGenomicFeature.indicator = TRUE,
                   alpha.set, SimilarityMatrix,
                   PCA.indicator = FALSE, PCA.feature, 
                   comp.num, fun = c("logcosh", "exp"), fun.alpha = 1,
                   maxit = 200, tol = 1e-04, verbose = FALSE, w.init = NULL) {
  
  library(Matrix)
  
  dd <- dim(InputData)
  d <- dd[dd != 1L]
  if (length(d) != 2L)
    stop("Input data must be matrix-conformal")
  InputData <- if (length(d) != length(dd)) matrix(InputData, d[1L], d[2L])
  else as.matrix(InputData)
  
  if (nrow(InputData) >= ncol(InputData)) {
    message("Transpose input data")
    InputData <- t(InputData)
  }
  
  if (fun.alpha < 1 || fun.alpha > 2)
    stop("The parameter fun.alpha must be in range [1,2]")
  fun <- match.arg(fun)
  
  
  if (PCA.indicator) {
    if (verbose) message("PCA analysis")
    if (missing(PCA.feature) || PCA.feature > min(nrow(InputData), ncol(InputData))) {
      stop("The features number in PCA analysis is incorrect")
    } else {
      use.data <- pca.cov(InputData, topNfeat = PCA.feature)
    }
  } else {
    use.data <- InputData
  }
  
  
  row.num <- nrow(use.data)
  col.num <- ncol(use.data)
  
  if(missing(comp.num)) {
    message("parameter 'comp.num' missing: set to ", min(row.num, col.num))
    comp.num <- min(row.num, col.num)
  }
  if (comp.num > min(row.num, col.num)) {
    message("'comp.num' is too large: reset to ", min(row.num, col.num))
    comp.num <- min(row.num, col.num)
  }
  
  if(is.null(w.init))
    w.init <- matrix(rnorm(comp.num ^ 2), comp.num, comp.num)
  else {
    if(!is.matrix(w.init) || length(w.init) != (comp.num ^ 2))
      stop("w.init is not a matrix or is the wrong size")
  }
  
  if (FusionGenomicFeature.indicator) {
    if (verbose) message("Remove the correlation between samples or conditions")
    if (verbose) message("Calculate Matrix inverse!")
    if(missing(alpha.set)) {
      message("parameter 'alpha.set' missing: set to ", 0.5)
      alpha.set <- 0.5
    }
    
    if(alpha.set == 1) {
      use.inverse <- as.matrix(diag(col.num))
    } else {
      unit.matrix <- diag(col.num)
      use.matrix <- alpha.set * unit.matrix + (1 - alpha.set) * SimilarityMatrix
      incidence.colsum <- apply(SimilarityMatrix, 2, sum)
      incidence.index <- which(incidence.colsum == 0)
      use.matrix[, incidence.index] <- unit.matrix[, incidence.index]
      
      rm(unit.matrix)
      rm(SimilarityMatrix)
      rm(incidence.colsum)
      rm(incidence.index)
      gc()
      
      use.matrix.lu <- expand(lu(as.matrix(use.matrix)))
      use.inverse <- solve(use.matrix.lu$U) %*% solve(use.matrix.lu$L) %*% t(use.matrix.lu$P)
      use.inverse <- apply(use.inverse, 2, as.numeric)
      
    }
    data.incidence <- use.data %*% use.inverse
    rm(use.inverse)
    gc()
  } else {
    data.incidence <- use.data
  }
  
  if (verbose) message("Centering")
  data.incidence <- t(scale(t(data.incidence), scale = FALSE))
  
  if (verbose) message("Whitening")
  cov.matrix <- data.incidence %*% t(data.incidence)/col.num
  cov.svd <- La.svd(cov.matrix)
  whiten.matrix <- diag(c(1 / sqrt(cov.svd$d))) %*% t(cov.svd$u)
  whiten.matrix <- matrix(whiten.matrix[1:comp.num, ], comp.num, row.num)
  whiten.data <- whiten.matrix %*% data.incidence
  
  # main loop of ICA analysis.
  if (verbose) message("parallel ICA analysis")
  result <- ica.par(whiten.data, tol = tol, fun = fun, fun.alpha = fun.alpha,
                    maxit = maxit, verbose = verbose, w.init = w.init)
  
  w <- result %*% whiten.matrix
  S <- w %*% data.incidence
  A <- t(w) %*% solve(w %*% t(w))
  if (PCA.indicator) {
    return(list(I = InputData, X = t(use.data), K = t(whiten.matrix), 
                W = t(result), A = t(A), S = t(S)))
  } else {
    return(list(X = InputData, K = whiten.matrix, 
                W = result, A = A, S = S))
  }
  if (verbose) message("Successful")
}


# PCA dimensionality reduction.
pca.cov <- function(input, topNfeat) {
  meanVals <- apply(input, 2, mean)
  DataAdjust <- sweep(input, 2, meanVals)
  covMat <- cov(DataAdjust)
  eigVects <- eigen(covMat)$vec
  redEigVects <- eigVects[, 1:topNfeat]
  lowDDataMat <- DataAdjust %*% redEigVects #将数据转换到低维新空间
  return(lowDDataMat)
}


# parallel Fast ICA.
ica.par <- function(X, tol, fun, fun.alpha, maxit, verbose, w.init) {
  Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
  p <- ncol(X)
  W <- w.init
  sW <- La.svd(W)
  #奇异值就是矩阵特征值的平方根！！！
  W <- sW$u %*% Diag(1 / sW$d) %*% t(sW$u) %*% W
  W1 <- W
  lim <- rep(1000, maxit)
  it <- 1
  if (fun == "logcosh") {
    if (verbose)
      message("Symmetric FastICA using logcosh approx. to neg-entropy function")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- tanh(fun.alpha * wx)
      v1 <- gwx %*% t(X)/p
      g.wx <- fun.alpha * (1 - (gwx) ^ 2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1 / sW1$d) %*% t(sW1$u) %*% W1
      
      # lim[it + 1] <- max(apply(abs(W1 - W), 2, sum))
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      
      W <- W1
      if (verbose)
        message("Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  if (fun == "exp") {
    if (verbose)
      message("Symmetric FastICA using exponential approx. to neg-entropy function")
    while (lim[it] > tol && it < maxit) {
      wx <- W %*% X
      gwx <- wx * exp( - (wx ^ 2) / 2)
      v1 <- gwx %*% t(X)/p
      g.wx <- (1 - wx ^ 2) * exp( - (wx ^ 2) / 2)
      v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
      W1 <- v1 - v2
      sW1 <- La.svd(W1)
      W1 <- sW1$u %*% Diag(1 / sW1$d) %*% t(sW1$u) %*% W1
      lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
      W <- W1
      if (verbose)
        message("Iteration ", it, " tol = ", format(lim[it + 1]))
      it <- it + 1
    }
  }
  W
}



#PCA.
pca.cov <- function(input, topNfeat){
  meanVals <- apply(input, 2, mean)
  DataAdjust <- sweep(input, 2, meanVals)
  covMat <- cov(DataAdjust)
  eigVects <- eigen(covMat)$vec
  redEigVects <- eigVects[, 1:topNfeat]
  lowDDataMat <- DataAdjust %*% redEigVects #将数据转换到低维新空间
  return(lowDDataMat)
}


GetPatterns <- function(ICA_result, pattern.var) {
  analysis.result <- as.matrix(t(ICA_result$S))
  up.index <- data.frame()
  bottom.index <- data.frame()
  up.var <- 1
  bottom.var <- 1
  for (variable in 1:ncol(analysis.result)) {
    temp <- analysis.result[, variable]
    temp.up <- which(temp >= pattern.var * var(temp))
    if (length(temp.up) != 0) {
      up.index[1:length(temp.up), up.var] <- temp.up
      up.var <- up.var + 1
    }
    
    temp.bottom <- which(temp <= - pattern.var * var(temp))
    if (length(temp.bottom) != 0) {
      bottom.index[1:length(temp.bottom), bottom.var] <- temp.bottom
      bottom.var <- bottom.var + 1
    }
  }
  up.index <- as.matrix(up.index)
  bottom.index <- as.matrix(bottom.index)
  return(list(up.index = up.index, bottom.index = bottom.index))
}