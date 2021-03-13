#' @title Class container of results obtained by spectra annotation
#'
#' @description
#'
#' `AnnotdSpectra` object provides a structure to the annotation results of a
#' query spectra  against a  reference spectra library. It must contain
#' reference spectra (with metadata) along with their reference compound
#' metadata, the query spectra (with metadata), and a hits table linking all
#' three items; the hits table also contains the distance metrics and other
#' variables that may define a hit (e.g. adduct assumed on the query spectrum)
#'
#' @details
#'
#' `AnnotdSpectra` objects should be created using the constructor function
#' `AnnotdSpectra` providing required spectra as a spectra object
#' (TODO:Reference).
#'
#' @slot hits cross-reference dataframe containing the hits along with their
#'   proposed adducts and common masses
#' @slot refCompound dataframe containing (reference) compound metadata
#' @slot qrySpectra spectra(n) with the query spectra
#' @slot refSpectra spectra(n) with reference spectra
#' @slot .properties inner info
#' @slot infoAnnotation variables used on the annotation process
#'
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
