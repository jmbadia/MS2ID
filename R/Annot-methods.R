#' @include Annot.R

#' @export
#' @describeIn Annot Getter for the refCompound slot
setMethod("refCompound", "Annot", function(x) x@refCompound)

#' @export
#' @describeIn Annot Getter for the qrySpectra slot
setMethod("qrySpectra", "Annot", function(x) x@qrySpectra)

#' @export
#' @describeIn Annot Getter for the refSpectra slot
setMethod("refSpectra", "Annot", function(x) x@refSpectra)

#' @export
#' @describeIn Annot Getter for the hits slot
setMethod("hits", "Annot", function(x) x@hits)

#' @export
#' @describeIn Annot Getter for the infoAnnotation slot
setMethod("infoAnnotation", "Annot", function(x) x@infoAnnotation)
