cal_we_score <- function(gene.vector) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  gene.vector <- as.vector(unique(na.omit(gene.vector)))
  go <- enrichGO(gene.vector, 'org.Hs.eg.db', ont = "ALL", pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID', readable = TRUE)
  if (is.null(go)) {
    we.score <- 0
  } else {
    if (nrow(go@result) == 0) {
      we.score <- 0
    } else {
      result <- matrix(nrow = nrow(go@result), ncol = 2)
      result[, 1] <- go@result[, 10]
      result[, 2] <- as.numeric(go@result[, 6])
      
      gene.get <- go@result[, 9]
      gene.get <- paste(gene.get[1:length(gene.get)], collapse = "/")
      gene.use.num <- length(unique(unlist(strsplit(gene.get, "[/]"))))
      
      result[, 1] <- result[, 1] / length(gene.vector)
      result[, 2] <- -log(result[, 2], 10)
      denominator <- (apply(result, 2, sum))[1]
      denominator <- denominator + ((length(gene.vector) - gene.use.num) / length(gene.vector))
      molecule <- sum(result[, 1] * result[, 2])
      we.score <- molecule / denominator
    }
  }
  return(we.score)
}