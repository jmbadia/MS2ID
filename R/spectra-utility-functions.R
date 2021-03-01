## Utility functions to preprocess query Spectra.

#' Remove invalid spectra
#'
#' '.validateSpectra' checks all spectra matrices and remove the invalid
#'  ones from the list of matrices and metadata. A spectra matrix is
#'  considered valid when it has two rows and at list one column.
#'
#' @param QRY object with the structure of .loadSpectra() return
#'
#' @return The QRY object with any invalid spectra removed, from the metadata
#' but also from the list of matrices
#' @noRd

.validateSpectra <- function(QRY){
  invalidMatrx <- vapply(QRY$Spectra, function(x) {
    if(class(x)[1] == "matrix") nrow(x) != 2 | ncol(x) < 1
    else TRUE
  }, FUN.VALUE = T)

  if(all(invalidMatrx)){
    stop("Query samples do not have valid spectra")
  } else if (any(invalidMatrx)){
    QRY <- .pruneSpectra(QRY, names(QRY$Spectra)[invalidMatrx])
    }
  return(QRY)
}


#' remove spectra from QRY
#'
#' '.pruneSpectra' removes spectra based on its id from the metadata
#' but also from the list of matrices
#'
#' @param QRY object with the structure of .loadSpectra() return
#' @param idSpectra2remove vector of spectrum identfiers (integers)
#'  pointing out the spectra to be removed
#'
#' @noRd

.pruneSpectra <- function(QRY, idSpectra2remove){
    QRY$Metadata <- QRY$Metadata[!(QRY$Metadata$idSpectra %in% idSpectra2remove)
                                 , ]

    keepSpctra <- !(names(QRY$Spectra) %in% idSpectra2remove)
    QRY$Spectra <- QRY$Spectra[keepSpctra]
    return(QRY)
}

#' Bin spectra
#'
#' '.binSpectra' rounds mz masses of spectra and merge those with the same
#' resulting value.
#'
#' @param spectraList list where every item is a spectrum. Every spectrum is a
#' matrix with two rows (named 'mass-charge' and 'intensity') and a column
#' for every mass.
#'
#' @param decimals2bin integer or numeric indicating the number of decimal places to be used on the round.
#'
#' @return The spectraList list binned
#' @noRd

.binSpectra <- function(spectraList, decimals2bin){
    #check argument types
    reqClasses <- c(decimals2bin="integer")
    .checkTypes(as.list(match.call(expand.dots=FALSE))[-1], reqClasses)

    if (decimals2bin < 0)
      stop("'decimals2bin' is expected to be a natural number")

    rmbrNames<-rownames(spectraList[[1]])
    # round spectral masses and sum their intensities up
    # if the resulting mass matches
    spectraList <- pbapply::pblapply(spectraList, function(x){
        x["mass-charge", ]<-round(x["mass-charge", ], decimals2bin)
        a <- t(x)
        x <- t(stats::aggregate(a[ ,"intensity"],
                                by=list(a[ ,"mass-charge"]), sum))
        rownames(x) <- rmbrNames
        return(x)
    })
    return(spectraList)
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
  spectraList[[as.character(idSpectrum)]][rowName, ]
}




