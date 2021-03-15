#' @include AnnotdSpectra.R

#' @describeIn AnnotdSpectra Getter for the refCompound slot
setMethod("refCompound", "AnnotdSpectra", function(x) x@refCompound)

#' @describeIn AnnotdSpectra Getter for the qrySpectra slot
setMethod("qrySpectra", "AnnotdSpectra", function(x) x@qrySpectra)

#' @describeIn AnnotdSpectra Getter for the refSpectra slot
setMethod("refSpectra", "AnnotdSpectra", function(x) x@refSpectra)

#' @describeIn AnnotdSpectra Getter for the hits slot
setMethod("hits", "AnnotdSpectra", function(x) x@hits)

#' @describeIn AnnotdSpectra Getter for the infoAnnotation slot
setMethod("infoAnnotation", "AnnotdSpectra", function(x) x@infoAnnotation)
