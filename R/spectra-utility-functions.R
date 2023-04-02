## Utility functions to preprocess Spectra.

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
    if(is(x, "matrix")) nrow(x) != 2 | ncol(x) < 1
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

#' Bin ONLY STATICALLY spectra
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
    reqClasses <- c(decimals2bin = "integer")
    .checkTypes(c(as.list(environment())), reqClasses)
    if (decimals2bin < 0)
      stop("'decimals2bin' is expected to be a natural number")
    lapply(spectraList, function(x) .binSpectrum(x, decimals2bin))
}

#' Bin spectrum statically or dynamically
#'
#' '.binSpectrum' bins the masses of a given spectrum either statically or
#' dynamically based on which argument is chosen, decimals2bin or massError. The
#' resulting intensity of the merged mz values is the sum of their respective
#' intensities. In the dynamic case, binSpectrum identifies each fragment and
#' searches for its corresponding fragments with similar masses (based on mass
#' error) in descending order of intensity.
#'
#' @param spectrum a two-row matrix (named 'mass-charge' and 'intensity') with a
#'   column for each mass value..
#'
#' @param decimals2bin An integer or numeric value indicating the number of
#'   decimal places to be used for rounding.
#'
#' @param massError An integer or numeric value indicating the ppm error
#'   considered for merging mz values.
#'
#' @return The spectrum binned
#' @noRd

.binSpectrum <- function(spectrum, decimals2bin, massError){
    if(ncol(spectrum) == 1) return(spectrum)

    if (missing(massError)) { #bin statically
        spectrum["mass-charge", ] <- round(spectrum["mass-charge", ],
                                           decimals2bin)
        a <- t(spectrum)
        spectrum <- t(stats::aggregate(a[ ,"intensity"],
                                       by = list(a[ ,"mass-charge"]), sum))
        rownames(spectrum) <- c("mass-charge", "intensity")
    } else{ # or bin dinamically
        #sort matrix according m/z
        spectrum <- spectrum[, order(spectrum['mass-charge',])]
        #obtain positions sorted by intensity
        colsByIntens <- sort(spectrum['intensity',], index.return = T,
                             decreasing = T)$ix # cols sorted by intensity
        cols2Keep <- colsByIntens
        for(col2kp in colsByIntens){ #ordered by Intens, we look for equal mz
            # if this mz is not been removed by a former loop
            if(col2kp %in% cols2Keep){
                mzSimilars <- MS2ID:::.posWhere(spectrum['mass-charge',],
                                                spectrum['mass-charge', col2kp],
                                                massError)
                if(length(mzSimilars) > 1){
                    #sum intensities
                    spectrum['intensity', col2kp] <- sum(spectrum['intensity',
                                                                  mzSimilars])
                    # put NA in similars mz
                    cols2Keep[cols2Keep %in% setdiff(mzSimilars, col2kp)] <- NA
                }
            }
        }
        spectrum <- spectrum[, sort(na.omit(cols2Keep)), drop=FALSE]
        #returns merged fragments sorted by mz
    }
    return(spectrum)
}

#' Extract a spectrum from a spectra list
#'
#' @param spectraList spectra list. spectra id are in names(list)
#' @param idSpectrum integer(1) or character(1) id of the spectrum to extract
#' @param rowName character(1). Optionally, 'mass-charge' or 'intensity'
#'  row can be selected
#'
#' @noRd
.getFragments <- function(spectraList, idSpectrum, rowName){
  if(length(idSpectrum)!=1) stop("Only a idSpectrum is supported")
  if(missing(rowName)){
    spectraList[[as.character(idSpectrum)]]
  }else{
    spectraList[[as.character(idSpectrum)]][rowName, ]
  }
}




