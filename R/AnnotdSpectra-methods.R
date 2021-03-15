#' @include AnnotdSpectra.R

#' @export
#' @describeIn AnnotdSpectra Getter for the refCompound slot
setMethod("refCompound", "AnnotdSpectra", function(x) x@refCompound)

#' @export
#' @describeIn AnnotdSpectra Getter for the qrySpectra slot
setMethod("qrySpectra", "AnnotdSpectra", function(x) x@qrySpectra)

#' @export
#' @describeIn AnnotdSpectra Getter for the refSpectra slot
setMethod("refSpectra", "AnnotdSpectra", function(x) x@refSpectra)

#' @export
#' @describeIn AnnotdSpectra Getter for the hits slot
setMethod("hits", "AnnotdSpectra", function(x) x@hits)

#' @export
#' @describeIn AnnotdSpectra Getter for the infoAnnotation slot
setMethod("infoAnnotation", "AnnotdSpectra", function(x) x@infoAnnotation)
