#' @title S4 class containing annotation results
#'
#' @aliases Annot-class
#'
#' @description
#'
#' `Annot` object structures the annotation results returned by the
#' \code{\link{annotate}} function. It contains reference spectra (with
#' metadata) along with their reference compound metadata, the query spectra
#' (with metadata), and a hits table linking all three items; the hits table
#' also contains the distance metrics and other Variables that may define a hit
#' (e.g. adduct assumed on the query spectrum).
#'
#' An \code{Annot} object is created internally by the \code{\link{annotate}}
#' function, but its content can be accessed externally by using: \itemize{\item
#' Getters of the object. See next section \item \code{\link{MS2IDgui}}
#' function: A GUI interface that browses visually the Annot content.\item
#' \code{\link{export2xlsx}} function: Export all the data to an xlsx file.}
#'
#' More info in the corresponding
#' \href{https://jmbadia.github.io/MS2ID/articles/annotate.html#annot-object}{vignette}.
#'
#' @section Getters: To obtain the content of any Annot's slot, please use the
#'   following methods \itemize{ \item hits: returns a cross-reference dataframe
#'   containing the hits along with their proposed adducts and common masses.
#'   \item qrySpectra(object): returns a \code{\link[Spectra]{Spectra}} object
#'   \insertCite{Spectra}{MS2ID} containing the query spectra with hits. \item
#'   refSpectra(object): returns a \code{\link[Spectra]{Spectra}} object with
#'   successful reference spectra. \item refCompound(object): returns a
#'   dataframe containing (reference) compound metadata of successful reference
#'   spectra. \item infoAnnotation(object): returns the variables used on the
#'   annotation process }
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#' @seealso \code{\link{annotate}} function, and post annotation tools
#'   \code{\link{MS2IDgui}} and \code{\link{export2xlsx}}.
#' @example man/examples/loadMS2ID.R
#' @example man/examples/selectQuerySpectra.R
#' @examples
#' ## ANNOTATION ---
#' library(MS2ID)
#' MS2IDlib <- MS2ID(MS2IDFolder)
#' annotResult <- annotate(QRYdata = queryFolder, MS2ID = MS2IDlib)
#'
#' ## ANNOTATION HANDLING---
#' #merge hits and compound info
#' result <- merge(x = hits(annotResult), y = refCompound(annotResult),
#'                by.x = "idREFcomp", by.y = "id", all.y = FALSE)
#' #Subset spectra and metadata considering first hit query spectra
#' library(Spectra)
#' idQRYspect_1 <- result$idQRYspect[1]
#' result_1 <- dplyr::filter(result, idQRYspect  == idQRYspect_1)
#' qrySpct_1 <- qrySpectra(annotResult)
#' qrySpct_1 <- qrySpct_1[qrySpct_1$id %in% result_1$idQRYspect]
#' refSpct_1 <- refSpectra(annotResult)
#' refSpct_1 <- refSpct_1[refSpct_1$id %in% result_1$idREFspect]
#' @importFrom Rdpack reprompt
#' @references \insertRef{Spectra}{MS2ID}
#' @name Annot
NULL

#' The Annot class
#'
#' An S4 class to represent annotation results
#'
#' @slot refCompound  data frame containing metadata of the reference compounds
#'   present in the \code{hits} table.
#' @slot qrySpectra \code{\link[Spectra]{Spectra}} object
#'   (\href{https://www.bioconductor.org/packages/release/bioc/html/Spectra.html}{Spectra
#'    package}) containing both successful query and consensus spectra (and
#'   their source spectra).
#' @slot refSpectra \code{\link[Spectra-class]{Spectra}} object
#'   (\href{https://www.bioconductor.org/packages/release/bioc/html/Spectra.html}{Spectra
#'    package}) with the reference spectra present in the \code{hits} table.
#' @slot hits cross-reference data frame containing the annotation hits, the id
#'   of the spectra and compounds and the columns \code{propAdduct} (the
#'   proposed adduct that would match the query precursor mass with the neutral
#'   mass of the reference compound) and  \code{cmnMasses} (number of fragments
#'   in common between the query and the reference spectrum)
#' @slot infoAnnotation variables used on the \code{annotate} function.
#' @slot .properties A list with properties of the class
#'
# #' @name Anot-class
# #' @docType class
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#'
#' @importFrom methods new
#' @importClassesFrom Spectra Spectra
#' @exportClass Annot
#' @noRd
.Annot <- setClass("Annot",
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
#' @noRd
setValidity("Annot", function(object) {
    slotN <- slotNames(object)
    nullV <- vapply(slotN, function(x) is.null(slot(object, x)), FUN.VALUE = T)
    if (!any(nullV))
        .validAnnot(object@dbcon, object@spectracon, object@mzIndexcon)
    else TRUE
})

#TODO: check structural integrity (required columns, ...)
.validAnnot <- function(q, rs, rf, h, i){
    txt <- character()
    if (length(txt)) txt else TRUE
}

#' Creator of the Annot class
#' @noRd
Annot <- function(qrySpectra, refSpectra, refCompound, hits, infoAnnotation) {
    #check if any missing arg
    argsDef <- ls()
    argsDeclr <- names(as.list(match.call())[-1])

    if (any(!(argsDef %in% argsDeclr)))
        stop(paste("Argument/s", paste(setdiff(argsDef, argsDeclr),
                                       collapse=", "),"is/are required"))

    res <- .validAnnot(q = qrySpectra, rs= refSpectra, rf = refCompound,
                               h = hits, i = infoAnnotation)
    if (is.character(res))
        stop(res)
    finalObj <- .Annot(qrySpectra = qrySpectra, refSpectra= refSpectra,
                               refCompound = refCompound, hits = hits,
                               infoAnnotation = infoAnnotation)
    return(finalObj)
}
