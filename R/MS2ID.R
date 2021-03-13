#' @title High-throughput MS/MS annotation with an unique in house database
#'
#' @description
#'
#' `MS2ID` package provides a class (MS2ID) to encapsulate a MS2 spectra
#' database along with its metadata. Its internal structure allows to annotate
#' query spectra - with the 'annotate' function - using big spectra databases at
#' high speed and low RAM requirements (typically 100 query spectra/min against
#' a 1M5 spectra library).
#' MS2ID class uses a SQLite database with the metadata
#' and bigmemory files to  store the spectra (i.e. peaks matrices) and their
#' mass-charge index.
#'
#' @details
#'
#' `MS2ID` objects should be created using the constructor function `MS2ID`
#' providing any file name (with path) related to the database (i.e. SQL file
#' name or any bigmemory file names)
#'
#' @section General functions:
#'
#'   - `MS2ID()`: Constructor of a MS2ID object
#'   - `annotate()`: function to annotate query spectra against a MS2ID library
#'
#' @author Josep M. Badia
#'
#' @seealso
#'
#' [annotate()] for the function to annotate query spectra against the MS2ID
#' database
#'
#' @docType package
#' @name MS2ID
NULL

#' S4 class containing a spectra library
#'
#' `MS2ID` objects provide a structure to encapsulate ready-to-use MS2 spectra
#' database along with its metadata. This internal structure allows handling
#' big spectra databases and annotating query spectra at high speed and low
#' RAM requirement.
#' The metadata are stored in a SQLite database; the spectra (i.e.
#' peaks matrices) and a mz index in bigmemory external files.
#'
#' @slot dbcon DBIConnection to metadata and pointers to big.matrix
#'   data
#' @slot spectracon a big.matrix object containing spectra
#' @slot mzIndexcon a big.matrix object containing a mass-charge index of the
#'   spectra
#' @slot .properties inner info
#'
#' @importFrom methods new
#'
#' @importClassesFrom bigmemory big.matrix
#'
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

#' MS2ID constructor
#'
#' @param x `character(1)` defining the directory path containing
#'   the MS2ID db files ('metadataDB.db' and bigmemory files)
#'
#' @export
#'
#' @importFrom DBI dbDriver
#' @importFrom RSQLite dbConnect
MS2ID <- function(dir) {
    if (missing(dir))
        stop("Argument 'dir' is required")
    if (!is.character(dir))
        stop("Argument 'dir' must be a character")
    if (!dir.exists(dir))
        stop(paste(basename(dir), "is not a valid directory"))
    dbFiles <- c("metadataDB.db", "mzIndex_body.bin", "mzIndex_body.desc",
                 "spectra_body.bin", "spectra_body.desc")
    if (!all(file.exists(file.path(dir, dbFiles))))
        stop(paste0(basename(dir),
                    " does not contain all the necessary files (",
                    paste(dbFiles, collapse=", "), ")"))

    SQLx <- dbConnect(dbDriver("SQLite"),
                      dbname = file.path(dir, "metadataDB.db"))
    bigM_mzI <- bigmemory::attach.big.matrix("mzIndex_body.desc",
                                             backingpath=dir)
    bigM_s <- bigmemory::attach.big.matrix("spectra_body.desc", backingpath=dir)
    res <- .validMS2ID(db=SQLx, mzI=bigM_mzI, spctr=bigM_s)
    if (is.character(res))
        stop(res)
    ms2idObj <- .MS2ID(dbcon=SQLx, spectracon= bigM_s, mzIndexcon= bigM_mzI)
    return(ms2idObj)
}
