## Utility functions to preprocess query Spectra.

#' Remove invalid spectra
#'
#' '.validateSpectra' checks all spectra matrices and remove the invalid
#'  ones from the list of matrices and metadata. A spectra matrix is
#'  considered valid when it has two rows and at list one column.
#'
#' @param DB object with the structure of .loadSpectra() return
#'
#' @return The db object with any invalid spectra removed, from the metadata
#' but also from the list of matrices
#' @noRd

.validateSpectra <- function(DB){
  invalidMatrx <- vapply(DB$Spectra$spectra, function(x) {
    if(class(x)[1] == "matrix") nrow(x) != 2 | ncol(x) < 1
    else TRUE
  }, FUN.VALUE = T)

  if(all(invalidMatrx)){
    stop("Database does not have valid spectra")
  } else if (any(invalidMatrx)){
    DB <- .pruneSpectra(DB, DB$Spectra$idSpectra[invalidMatrx])
    }
  return(DB)
}


#' remove spectra from DB
#'
#' '.pruneSpectra' removes spectra based on its id from the metadata
#' but also from the list of matrices
#'
#' @param DB object with the structure of .loadSpectra() return
#' @param idSpectra2remove vector of spectrum identfiers (integers)
#'  pointing out the spectra to be removed
#'
#' @noRd

.pruneSpectra <- function(DB, idSpectra2remove){
    DB$Metadata <- DB$Metadata[!DB$Metadata$idSpectra %in% idSpectra2remove,]

    posInvSpctra <- which(DB$Spectra$idSpectra %in% idSpectra2remove)
    DB$Spectra$idSpectra <- DB$Spectra$idSpectra[-posInvSpctra]
    DB$Spectra$spectra <- DB$Spectra$spectra[-posInvSpctra]
    return(DB)
}

#' Bin spectra
#'
#' '.binSpectra' rounds mz masses of spectra and merge those with the same
#' resulting value.
#'
#' @param spectra list where every item is a spectrum. Every spectrum is a
#' matrix with two rows (named 'mass-charge' and 'intensity') and a column
#' for every mass.
#'
#' @param decimals2bin integer or numeric indicating the number of decimal places to be used on the round.
#'
#' @return The spectra list binned
#' @noRd

.binSpectra <- function(spectra, decimals2bin){
    #check argument types
    reqClasses <- c(decimals2bin="integer")
    .checkTypes(as.list(match.call(expand.dots=FALSE))[-1], reqClasses)

    if (decimals2bin < 0)
      stop("'decimals2bin' is expected to be a natural number")

    rmbrNames<-rownames(spectra[[1]])
    # round spectral masses and sum their intensities up
    # if the resulting mass matches
    spectra <- pbapply::pblapply(spectra, function(x){
        x["mass-charge", ]<-round(x["mass-charge", ], decimals2bin)
        a <- t(x)
        x <- t(stats::aggregate(a[ ,"intensity"],
                                by=list(a[ ,"mass-charge"]), sum))
        rownames(x) <- rmbrNames
        return(x)
    })
    return(spectra)
}


#' Extract a spectrum from a spectra list
#'
#' @param spectraList spectra list. spectra id are in names(list)
#' @param idSpectrum integer(1) or character(1) id of the spectrum to extract
#' @param rowName character(1). Optionally, 'mass-charge' or 'intensity'
#'  row can be selected
#'
#' @return
#' @noRd
.getFragments <- function(spectraList, idSpectrum, rowName=""){
  if(length(idSpectrum)!=1) stop("Only a idSpectrum is supported")
  spectraList[[as.character(idSpectrum)]][rowName,]
}




