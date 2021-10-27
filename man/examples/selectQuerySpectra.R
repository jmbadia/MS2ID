
## SELECT QUERY SPECTRA ---
## Decompress the query mzML files that come with MS2ID
queryFile <- system.file("extdata/QRYspectra.zip",
                         package = "MS2ID",
                         mustWork = TRUE)
queryFolder <- file.path(tempdir(), "QRYspectra")
library(utils)
unzip(queryFile, exdir = queryFolder)
