## 1. Obtain a MS2ID library
library(CompoundDb)
wrkDir <- tempdir()
MoNAsubset <- system.file("extdata/MoNAsubset.sdf.gz",
                          package = "MS2ID")
cmps <- suppressWarnings(compound_tbl_sdf(MoNAsubset))
spctr <- msms_spectra_mona(MoNAsubset, collapsed = TRUE)
spctr$predicted <- FALSE
metad <- data.frame(
    name = c("source", "url",
             "source_version","source_date"),
    value = c("MoNA",
              "https://mona.fiehnlab.ucdavis.edu/downloads",
              "v1", '07-09')
)
cmpdbDir <- file.path(wrkDir, "cmpdbDir")
if(dir.exists(cmpdbDir)){
    do.call(file.remove,
            list(list.files(cmpdbDir, full.names = TRUE)))
}else{
    dir.create(cmpdbDir)
}
db_file <- createCompDb(cmps, metadata = metad,
                        msms_spectra = spctr,
                        path = cmpdbDir)
cmpdb <- CompDb(db_file)
library(MS2ID)
MS2IDdirectory <- createMS2ID(cmpdb = cmpdb,
                              overwrite = TRUE, path = wrkDir)
MS2IDlib <- MS2ID(MS2IDdirectory)

## 2. Select query spectra
queryFile <- system.file("extdata/Met1.zip",
                         package = "MS2ID")
queryFolder <- file.path(wrkDir, "QRYspectra")
library(utils)
unzip(queryFile, exdir = queryFolder)

## 3.1 Simple annotate
annotResult <- annotate(QRYdata = queryFolder,
                        MS2ID = MS2IDlib)
## 3.2 Add external function to annotate
foo <- function(finalMatrix){
  vector1 <- finalMatrix[1,]
  vector2 <- finalMatrix[2,]
  CosplusOne <- 1+ suppressMessages(
    philentropy::distance(rbind(vector1, vector2),
                          method = "cosine")
  )
  names(CosplusOne) <- "CosplusOne"
  return(CosplusOne)
}
annotate(QRYdata = queryFolder, MS2ID = MS2IDlib,
          metrics = c("cosine", "fidelity", "topsoe"),
          metricsThresh = c(0.8, 0.6, 0.6),
          metricFUN = foo, metricFUNThresh = 1.8)
