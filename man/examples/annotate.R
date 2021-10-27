
## ANNOTATE ---
## 3.1 Simple annotation
library(MS2ID)
MS2IDlib <- MS2ID(MS2IDFolder)
annotResult <- annotate(QRYdata = queryFolder,
                        MS2ID = MS2IDlib)
## 3.2 Annotation with different metrics
## external function to annotate
foo <- function(finalMatrix){
  vector1 <- finalMatrix[1, ]
  vector2 <- finalMatrix[2, ]
  CosplusOne <- 1 + suppressMessages(
    philentropy::distance(rbind(vector1, vector2),
                          method = "cosine")
  )
  names(CosplusOne) <- "CosplusOne"
  return(CosplusOne)
}
annotResult <- annotate(
  QRYdata = queryFolder, MS2ID = MS2IDlib,
  metrics = c("cosine", "fidelity", "topsoe"),
  metricsThresh = c(0.8, 0.6, 0.6),
  metricFUN = foo, metricFUNThresh = 1.8)
