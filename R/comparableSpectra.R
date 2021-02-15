#' Obtain reference spectra comparable to a query spectrum
#' 
#' This function identifies reference spectra that have peaks in common with a 
#' query spectrum by consulting a mz index. Arguments allow modulate 
#' requirements of minimum number of common peaks (cmnPeaks) or their 
#' relative weigth (cmnTopPeaks).
#'
#' @param QRYmz numeric vector with mz values, presumably from a query spectrum.
#' @param cmnPeaks Filter by fragment argument. Resulting REF spectra must
#'  have at least a number N (cmnPeaks value) of common fragments
#'   with UNK spectra.
#' @param cmnTopPeaks REF library subsetting argument. Only REF spectra with
#'  spectra whose top N (topMassesNum value) fragments share at list one mz 
#'  with UNK fragments.
#'
#' @return A numeric vector with the id number of the reference spectra that 
#' meet the requirements.
.comparableSpectra <- function(QRYmz, cmnPeaks=2, cmnTopPeaks=5){
  #positions of lstmzIdx[[n]] with mz match
  pos_mzIndx <- lapply(DB$mzIdx$lstmzIdx, function(x) x %in% QRYmz)
  #In case no REF spectra has such QRYmz, EXIT
  if(!any(unlist(pos_mzIndx))) return(NULL)
  ## obtain spectra that contains any QRYmz
  resultSpectra <- lapply(seq_along(DB$mzIdx$lstmzIdx), function(topPos) 
    unlist(DB$mzIdx$spectrIdx[[topPos]][pos_mzIndx[[topPos]]]))
  ## take spectra with almost "cmnPeaks" number of QRYmz
  with_cmnPeaks <- which(tabulate(unlist(resultSpectra)) >= cmnPeaks)
  # intersect with spectra where at list one mz is among the top n peaks
  resultSpectra <- with_cmnPeaks[with_cmnPeaks %in% 
                                   unique(unlist(resultSpectra
                                                 [seq_len(cmnTopPeaks)]))]
  return(resultSpectra)
}