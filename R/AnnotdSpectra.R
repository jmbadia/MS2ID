#' @name AnnotdSpectra
#'
#' @title S4 class containing annotation results
#'
#' @aliases AnnotdSpectra-class
#'
#' @description
#'
#' `AnnotdSpectra` object structures the annotation results. It must contain
#' reference spectra (with metadata) along with their reference compound
#' metadata, the query spectra (with metadata), and a hits table linking all
#' three items; the hits table also contains the distance metrics and other
#' variables that may define a hit (e.g. adduct assumed on the query spectrum)
#'
#' @section Getters: To obtain the content of any AnnotdSpectra's slot, please
#'   use the following methods \itemize{ \item hits(object): returns a
#'   cross-reference dataframe containing the hits along with their proposed
#'   adducts and common masses. \item qrySpectra(object): returns a
#'   \code{\link[Spectra]{Spectra}} object (see
#'   \href{https://www.bioconductor.org/packages/release/bioc/html/Spectra.html}{Spectra
#'    package}) containing the query spectra with hits. \item
#'   refSpectra(object): returns a \code{\link[Spectra]{Spectra}} object with
#'   successful reference spectra. \item refCompound(object): returns a
#'   dataframe containing (reference) compound metadata of successful reference
#'   spectra. \item infoAnnotation(object): returns the variables used on the
#'   annotation process }
NULL

#' @importFrom methods new
#' @importClassesFrom Spectra Spectra
#' @exportClass AnnotdSpectra
.AnnotdSpectra <- setClass("AnnotdSpectra",
                    slots = c(refCompound = "data.frame",
                              qrySpectra = "Spectra",
                              refSpectra = "Spectra",
                              hits = "data.frame",
                              infoAnnotation = "list",
                              .properties = "list"),
                    prototype = list(refCompound = NULL,
                                     qrySpectra = NULL,
                                     refSpectra = NULL,
                                     hits = NULL,
                                     infoAnnotation = NULL,
                                     .properties = list()
                                     ))

#' @importFrom methods validObject
setValidity("AnnotdSpectra", function(object) {
    slotN <- slotNames(object)
    nullV <- vapply(slotN, function(x) is.null(slot(object, x)), FUN.VALUE = T)
    if (!any(nullV))
        .validAnnotdSpectra(object@dbcon, object@spectracon, object@mzIndexcon)
    else TRUE
})

#TODO: check structural integrity (required coliumns, ...)
.validAnnotdSpectra <- function(q, rs, rf, h, i){
    txt <- character()
    if (length(txt)) txt else TRUE
}


AnnotdSpectra <- function(qrySpectra, refSpectra, refCompound, hits,
                              infoAnnotation) {
    #check if any missing arg
    argsDef <- ls()
    argsDeclr <- names(as.list(match.call())[-1])

    if (any(!(argsDef %in% argsDeclr)))
        stop(paste("Argument/s", paste(setdiff(argsDef, argsDeclr),
                                       collapse=", "),"is/are required"))

    res <- .validAnnotdSpectra(q = qrySpectra, rs= refSpectra, rf = refCompound,
                               h = hits, i = infoAnnotation)
    if (is.character(res))
        stop(res)
    finalObj <- .AnnotdSpectra(qrySpectra = qrySpectra, refSpectra= refSpectra,
                               refCompound = refCompound, hits = hits,
                               infoAnnotation = infoAnnotation)
    return(finalObj)
}
