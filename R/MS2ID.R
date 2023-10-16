#' @title MS2ID class for high-throughput MS/MS annotation
#'
#' @description
#'
#'The `MS2ID` class encapsulates, a MS2 spectra database along with its
#'metadata. Its internal structure allows to annotate query spectra - with the
#'\code{\link{annotate}} function - using big spectra databases at high speed
#'and low RAM requirements (typically 100 query spectra/min against a 1M5
#'spectra library). MS2ID class uses a SQLite database with the metadata and
#'bigmemory files to store the spectra (i.e. peaks matrices) and their
#'mass-charge index. See
#'\href{https://jmbadia.github.io/MS2ID/articles/MS2ID.html}{vignette}.
#'
#'@details The \code{createMS2ID} function creates a MS2ID backend that is
#'  subsequently used by the \code{MS2ID} constructor, which creates the
#'  \code{MS2ID} object.
#'
#'@section General functions: \itemize{  \item createMS2ID(): Function to create
#'  MS2ID backends. \item MS2ID(): Constructor function for MS2ID objects.}
#'
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#' @seealso \code{\link{createMS2ID}} function.
#' @example man/examples/MS2IDcreation.R
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{CompoundDB}{MS2ID}
#' @aliases MS2ID-class
#' @name MS2ID
NULL

#' The MS2ID class
#'
#' An S4 class to represent a reference library.
#'
#' @slot dbcon A DBIConnection object with the MSQL data
#' @slot spectracon A big.matrix object with the spectra data
#' @slot mzIndexcon A big.matrix object with index of fragments
#' @slot .properties A list with properties of the class
#'
#' @name MS2ID-class
#' @docType class
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#'
#' @importFrom methods new
#' @importClassesFrom bigmemory big.matrix
#' @exportClass MS2ID
#' @noRd
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

#' Creator of the MS2ID class
#'
#' @param ms2idFolder character(1) with the directory's path of the MS2ID
#'   backend
#'
#' @rdname MS2ID
#' @importFrom DBI dbDriver
#' @importFrom RSQLite dbConnect dbSendQuery dbClearResult
#' @importFrom bigmemory attach.big.matrix
#' @export
MS2ID <- function(ms2idFolder) {
    if (missing(ms2idFolder))
        stop("Argument 'ms2idFolder' is required")
    if (!file.exists(ms2idFolder))
        stop(paste(ms2idFolder, "is not a valid directory"))

    dbFiles <- c("MS2ID.db", "mzIndex_body.bin", "mzIndex_body.desc",
                 "spectra_body.bin", "spectra_body.desc")
    if (!all(file.exists(file.path(ms2idFolder, dbFiles))))
        stop(glue::glue("
        {dirname(ms2idFolder)} does not contain all the necessary files: \\
        {glue::glue_collapse(dbFiles, ', ', last = ' and ')}
                        "))

    SQLx <- RSQLite::dbConnect(DBI::dbDriver("SQLite"),
                      dbname = file.path(ms2idFolder, "MS2ID.db"))
    RSQLite::dbClearResult(RSQLite::dbSendQuery(SQLx, "PRAGMA busy_timeout=5000;"));

    bigM_mzI <- bigmemory::attach.big.matrix("mzIndex_body.desc",
                                             backingpath = ms2idFolder)
    bigM_s <- bigmemory::attach.big.matrix("spectra_body.desc",
                                           backingpath = ms2idFolder)
    res <- .validMS2ID(db=SQLx, mzI=bigM_mzI, spctr=bigM_s)
    if (is.character(res))
        stop(res)
    ms2idObj <- new("MS2ID", dbcon=SQLx, spectracon= bigM_s, mzIndexcon= bigM_mzI)
    return(ms2idObj)
}

#' Subset method for Annot class
#'
#' @param object character: signature supported
#' @rdname MS2ID
#'
#' @importFrom DBI dbGetQuery
#' @importMethodsFrom methods show
#'
#' @exportMethod show
setMethod("show", "MS2ID",
          function(object) {
              crossRef <- dbGetQuery(object@dbcon,
                                     "SELECT * FROM crossRef_SpectrComp")
              cat("MS2ID spectra library (", class(object)[1L],
                  ") \nContains: ", length(unique(crossRef$ID_compound)),                        " compounds / ", length(unique(crossRef$ID_spectra)),
                  " spectra\n", sep = "")
              props <- dbGetQuery(object@dbcon,
                                  "SELECT * FROM dbOriginal")
              rdUpdate <- strptime(props$lastRawDataUpdate, "%Y%m%d_%H%M%S")
              rdUpdate <- as.character(rdUpdate)
              lastMod <- as.POSIXct(props$lastModification*24*60*60,
                                    origin = "1970-01-01", tz="UTC")
              lastMod <- as.character(lastMod)
              cat("Database original: ", props$ID_db,
                  "\nLast raw data update: ", rdUpdate,
                  "\nLast modification: ", lastMod, "\n", sep = "")
          })
