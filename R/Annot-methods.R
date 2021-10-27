#' @include Annot.R

#' Subset method for Annot class
#'
#' @param x character: signature supported
#' @rdname Annot
#' @exportMethod refCompound
setMethod("refCompound", "Annot", function(x) x@refCompound)

#' Subset method for Annot class
#'
#' @param x character: signature supported
#' @rdname Annot
#' @exportMethod qrySpectra
setMethod("qrySpectra", "Annot", function(x) x@qrySpectra)

#' Subset method for Annot class
#'
#' @param x character: signature supported
#' @rdname Annot
#' @exportMethod refSpectra
setMethod("refSpectra", "Annot", function(x) x@refSpectra)

#' Subset method for Annot class
#'
#' @param x character: signature supported
#' @rdname Annot
#' @exportMethod hits
setMethod("hits", "Annot", function(x) x@hits)

#' Subset method for Annot class
#'
#' @param x character: signature supported
#' @rdname Annot
#' @exportMethod infoAnnotation
setMethod("infoAnnotation", "Annot", function(x) x@infoAnnotation)

#' Subset method for Annot class
#'
#' @param object character: signature supported
#' @rdname Annot
#' @importMethodsFrom methods show
#' @exportMethod show
setMethod("show", "Annot",
          function(object) {
              hits <- hits(object)
              cat("Annot object with annotation results: \nContains ",
                  length(unique(hits$idQRYspect)),
                  " different query spectra annotated with ",
                  length(unique(hits$idREFcomp)),
                  " different compounds\n", sep = "")
              info <- object@infoAnnotation
              cat("Database original: ", info$MS2ID,
                  "\nQuery spectra source: ", info$QRYdata,
                  "\nAnnotation date: ", info$annotationTime,
                  "\nUse hits() for more parameters used in the annotation.",
                  "\n", sep = "")
          })
