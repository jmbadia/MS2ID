## OBTAIN MS2ID LIBRARY ---
## Decompress the MoNA subset that comes with MS2ID
MoNAsubset <- system.file("extdata/MoNAsubset.sdf.gz",
                          package = "MS2ID")
## Use CompoundDB to parse MoNA
library(CompoundDb)
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
cmpdbDir <- file.path(tempdir(), "cmpdbDir")
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
## Create the MS2ID backend
library(MS2ID)
MS2IDdirectory <- createMS2ID(cmpdb = cmpdb,
                              overwrite = TRUE, path = tempdir())
## Obtain the MS2ID object
MS2IDlib <- MS2ID(MS2IDdirectory)
