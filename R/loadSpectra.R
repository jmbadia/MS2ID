#' Load query spectra
#'
#' @param dirPath `character` defining the location of the mzml files.
#' @param msLevel `integer` or `numeric`, but a natural number. It subsets the
#'  spectra to be loaded to those with such msLevel. By default, msLevel=2L
#' @param nsamples `integer` or `numeric`, but a natural number. To speed up
#'  the identification process, the user can limit spectra to load to
#'  'nsamples'. The spectra are chosen at random.
#'
#' @return a list with 2 items, 'Metadata' and 'Spectra'. The former is a data frame with spectrum metadata. The latter is a list with two items, a list of spectra (under matrix form) and 'idSpectra' (vector of its spectra id.) Both 'Metadata' and 'Spectra' are linked using the 'idSpectra' variable.
#' @noRd

.loadSpectra <- function(dirPath, msLevel=2L, nsamples){
    #check types
    reqClasses <- c(dirPath="character", msLevel="integer", nsamples="integer")
    .checkTypes(as.list(match.call(expand.dots=FALSE))[-1], reqClasses)

    if (msLevel < 2)
        stop("'msLevel' is expected to be a natural number > 1")
    if (!missing(nsamples))
        if (nsamples < 1)
            stop("'nsamples' is expected to be a natural number > 0")

    #load files
    mzml_files <- dir(path = dirPath, pattern="*\\.mzml$",
                    ignore.case=TRUE, full.names = FALSE)
    if(length(mzml_files) < 1)
        stop("No mzml files found in 'dirPath' path")

  # Take and merge metadata & spectralMatrix from MS2 spectra
  mtdata <- data.frame()
  spctra <- list()

  #Capture spectra info from all mzmlfiles
  for(id_file in seq_along(mzml_files)){
    mzMLfile <- mzR::openMSfile(file.path(dirPath, mzml_files[id_file]))
    temp <- mzR::header(mzMLfile)
    #positions to catch
    pos2Catch <- temp$msLevel==msLevel
    #not an scan to read, jump to next file
    if(!any(pos2Catch)) next

    #load metadata & filename. Append filename to spectra metadata
    temp <- cbind(file= mzml_files[id_file], temp[pos2Catch,])

    # load spectra matrix
    temp_spctra <- lapply(which(pos2Catch), function(x) {
      m <- t(mzR::peaks(mzMLfile,x))
      rownames(m) <- c("mass-charge","intensity")
      return(m)
    })

    #merge info
    spctra <- c(temp_spctra, spctra)
    mtdata <- rbind(temp,mtdata)
    mzR::close(mzMLfile)
  }

#rslt$QRY_spectra <- QRY$Spectra$spectra[relevantQRYSpectra]
#names(rslt$QRY_spectra) <- QRY$Spectra$idSpectra[relevantQRYSpectra]

  # Apply same ID to metadata spectra and metadata matrix
  mtdata <- cbind(idSpectra=seq_len(nrow(mtdata)), mtdata)
  names(spctra) <- seq_len(nrow(mtdata))

  #spctra <- list(idSpectra=seq_len(nrow(mtdata)), spectra=spctra)

  #subset to a number 'nsamples' of random samples.
  # In this part of the code in order to maintain id traceability
  if(!missing(nsamples)){
    if(nsamples > nrow(mtdata))
        message(paste("'nsamples' is larger than the number of spectra",
                      "available. All spectra will be loaded."))
    else{
      #set.seed(1)
        spectra2subset <- sample(mtdata$idSpectra, nsamples)
        mtdata <- mtdata[mtdata$idSpectra %in% spectra2subset, ]
        spctra <- spctra[names(spctra) %in% spectra2subset]
      }
  }

  result <- list(Metadata=mtdata, Spectra=spctra)
  return(result)
}
