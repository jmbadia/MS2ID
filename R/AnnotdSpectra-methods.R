# Generics
setGeneric("refCompound", function(x) standardGeneric("refCompound"))
setGeneric("qrySpectra", function(x) standardGeneric("qrySpectra"))
setGeneric("refSpectra", function(x) standardGeneric("refSpectra"))
setGeneric("hits", function(x) standardGeneric("hits"))
setGeneric("infoAnnotation", function(x) standardGeneric("infoAnnotation"))

#accesors
setMethod("refCompound", "AnnotdSpectra", function(x) x@refCompound)
setMethod("qrySpectra", "AnnotdSpectra", function(x) x@qrySpectra)
setMethod("refSpectra", "AnnotdSpectra", function(x) x@refSpectra)
setMethod("hits", "AnnotdSpectra", function(x) x@hits)
setMethod("infoAnnotation", "AnnotdSpectra", function(x) x@infoAnnotation)
