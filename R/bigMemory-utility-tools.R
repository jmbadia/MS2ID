## Utility functions to deal with big memory

#' Obtain spectra from bigmemory
#'
#' @param MS2ID MS2ID object containing the reference spectra
#' @param idSpectra vector(n) with the id of the spectra to get
#'
#' @return A list. 'spectra' item contains all the spectra collapsed.
#' 'spectraPTR' item is a dataframe with id, numItems and startPos data;
#' Thaht data are need it to extract every spectrum individually by using
#' the .extractSpectrum() function
#' @noRd

.bufferSpectra <- function(MS2ID, idSpectra){
    SQLwhere <- .appendSQLwhere("id", idSpectra, mode="IN")
    spectraPTR <- .getSQLrecords(MS2ID,"*", "spectraPTR", SQLwhere)

    #import ALL spectra (on disk -> RAM) to avoid concurrent access to disk
    readPos <- unlist(lapply(seq_len(nrow(spectraPTR)), function(x)
        seq_len(spectraPTR$numItems[x]) + spectraPTR$startPos[x]))
    spectra <- MS2ID@spectracon[, readPos, drop=FALSE]

    #recalculate startPos considering spectra is a subset
    spectraPTR$startPos <- c(0, cumsum(head(spectraPTR$numItems,-1)))
    return(list(spectra=spectra, ptr=spectraPTR))
}

#' Extract a reference spectrum
#'
#' @param bufferedSpectra list containing buffered spectra; results from the
#' .bufferSpectra() function.
#' @param position integer defining the spectrum to extract
#'
#' @return A reference spectrum
#' @noRd

.getSpectrum <- function(bufferedSpectra, position){
    bufferedSpectra$spectra[, seq_along(bufferedSpectra$ptr$numItems[position])+
                                bufferedSpectra$ptr$startPos[position],
                       drop=FALSE]
}

