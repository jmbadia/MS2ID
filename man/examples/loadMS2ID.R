
## LOAD MS2ID LIBRARY ---
## Decompress the MS2ID library that comes with MS2ID
MS2IDzipFile <- system.file("extdata/MS2IDLibrary.zip",
                            package = "MS2ID",
                            mustWork = TRUE)
library(utils)
MS2IDFolder <- dirname(unzip(MS2IDzipFile,
                             exdir = tempdir()))[1]
