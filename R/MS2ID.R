#'@name MS2ID
#'
#'@title In-house database for high-throughput MS/MS annotation
#'
#'@aliases MS2ID-class
#'
#'@description
#'
#'The `MS2ID` class encapsulates a MS2 spectra database along with its metadata.
#'Its internal structure allows to annotate query spectra - with the
#'\code{\link{annotate}} function - using big spectra databases at high speed
#'and low RAM requirements (typically 100 query spectra/min against a 1M5
#'spectra library). MS2ID class uses a SQLite database with the metadata and
#'bigmemory files to  store the spectra (i.e. peaks matrices) and their
#'mass-charge index.
#'
#'@details
#'
#'`MS2ID` objects should be created using the constructor function `MS2ID`
#'providing the path where the library files are (i.e. SQL file
#'name or any bigmemory file names)
#'
#'@section General functions: \itemize{ \item `MS2ID()`: Constructor of a MS2ID
#'  object. \item `annotate()`: function to annotate query spectra against a
#'  MS2ID library }
#'
#'@author Josep M. Badia
NULL

#' @importFrom methods new
#'
#' @importClassesFrom bigmemory big.matrix
#' @exportClass MS2ID
.MS2ID <- setClass("MS2ID",
                    slots = c(dbcon = "DBIConnection",
                              spectracon= "big.matrix",
                              mzIndexcon= "big.matrix",
                              .properties = "list"),
                    prototype = list(dbcon = NULL,
                                     spectracon= NULL,
                                     mzIndexcon= NULL,
                                     .properties = list()
                                     ))

#' @importFrom methods validObject
setValidity("MS2ID", function(object) {
    if (!any(is.null(object@dbcon), is.null(object@spectracon),
             is.null(object@mzIndexcon)))
        .validMS2ID(object@dbcon, object@spectracon, object@mzIndexcon)
    else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery
.validMS2ID <- function(db, mzI, spctr) {
    txt <- character()
    tables <- dbListTables(db)
    required_tables <- c("dbOriginal", "crossRef_SpectrComp", "metaSpectrum",
                         "metaCompound","spectraPTR", "mzIndexPTR1",
                         "mzIndexPTR2","mzIndexPTR3", "mzIndexPTR4",
                         "mzIndexPTR5","mzIndexPTR6")
    got <- required_tables %in% tables
    if (!all(got))
        txt <- c(txt, paste0("Required tables ", paste0(required_tables[!got]),
                             "not found in the database"))
    # TODO: develop deeper validation
    # check integrity mysql db <-> bigM (through indexes and ncols of bigM)
    ## Check table columns
    if (length(txt)) txt else TRUE
}

#' @export
#' @rdname MS2ID
#' @importFrom DBI dbDriver
#' @importFrom RSQLite dbConnect
MS2ID <- function(ms2idFolder) { #must be the directory's path containing the MS2ID files
    if (missing(ms2idFolder))
        stop("Argument 'ms2idFolder' is required")
    if (!file.exists(ms2idFolder))
        stop(paste(ms2idFolder, "is not a valid directory"))

    ms2idFolder <- dirname(ms2idFolder)
    dbFiles <- c("MS2ID.db", "mzIndex_body.bin", "mzIndex_body.desc",
                 "spectra_body.bin", "spectra_body.desc")
    if (!all(file.exists(file.path(ms2idFolder, dbFiles))))
        stop(paste0(ms2idFolder,
                    " does not contain all the necessary files (",
                    paste(dbFiles, collapse=", "), ")"))

    SQLx <- dbConnect(dbDriver("SQLite"),
                      dbname = file.path(ms2idFolder, "MS2ID.db"))
    bigM_mzI <- bigmemory::attach.big.matrix("mzIndex_body.desc",
                                             backingpath = ms2idFolder)
    bigM_s <- bigmemory::attach.big.matrix("spectra_body.desc",
                                           backingpath = ms2idFolder)
    res <- .validMS2ID(db=SQLx, mzI=bigM_mzI, spctr=bigM_s)
    if (is.character(res))
        stop(res)
    ms2idObj <- .MS2ID(dbcon=SQLx, spectracon= bigM_s, mzIndexcon= bigM_mzI)
    return(ms2idObj)
}
