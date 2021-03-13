#' @include AnnotdSpectra.R

#' @export
setMethod("refCompound", "AnnotdSpectra", function(x) x@refCompound)

#' @export
setMethod("qrySpectra", "AnnotdSpectra", function(x) x@qrySpectra)

#' @export
setMethod("refSpectra", "AnnotdSpectra", function(x) x@refSpectra)

#' @export
setMethod("hits", "AnnotdSpectra", function(x) x@hits)

#' @export
setMethod("infoAnnotation", "AnnotdSpectra", function(x) x@infoAnnotation)
