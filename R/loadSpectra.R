#' Load query spectra
#'
#' @param mzmlData `character` defining either the directory containing the mzML files, or the files themselves.
#' @param msLevel `integer` or `numeric`, but a natural number. It subsets the
#'  spectra to be loaded to those with such msLevel. By default, msLevel=2L
#' @param nsamples `integer` or `numeric`, but a natural number. To speed up
#'  the identification process, the user can limit spectra to load to
#'  'nsamples'. The spectra are chosen at random.
#'
#' @param acquisitionNum vector of`integer` or `numeric`. The user can limit spectra based on their acquisitionNum
#'
#' @return a list with 2 items, 'Metadata' and 'Spectra'. The former is a data frame with spectrum metadata. The latter is a list with two items, a list of spectra (under matrix form) and 'idSpectra' (vector of its spectra id.) Both 'Metadata' and 'Spectra' are linked using the 'idSpectra' variable.
#' @noRd

.loadSpectra <- function(mzmlData, msLevel = 2L, nsamples, acquisitionNum){
    #check types
    reqClasses <- c(mzmlData="character", msLevel="integer", nsamples="integer", acquisitionNum = "integer")
    .checkTypes(as.list(match.call(expand.dots=FALSE))[-1], reqClasses)
    if (msLevel < 2)
        stop("'msLevel' is expected to be a natural number > 1")
    if (!missing(nsamples))
        if (nsamples < 1)
            stop("'nsamples' is expected to be a natural number > 0")
    if(identical(file_test("-d", mzmlData), T)){# one char and is dir
      #load files
      mzml_files <- dir(path = mzmlData, pattern="*\\.mzml$",
                        ignore.case=TRUE, full.names = T)
    }else{ #all char must be .mzml chars
      ismzml <- grepl("\\.mzml$", mzmlData, ignore.case = T)
      if(!all(ismzml))
        warning(paste("The following non mzML file/s will be ignored:",
                      paste(mzmlData[!ismzml], collapse=", ")))
      mzml_files <- mzmlData[ismzml]
    }
    if(length(mzml_files) < 1)
        stop("No valid mzML files found")

  # Take and merge metadata & spectralMatrix from MS2 spectra
  mtdata <- data.frame()
  spctra <- list()

  #Capture spectra info from all mzmlfiles #TODO: use lapply and do.call
  for(file in mzml_files){
    mzRobj <- mzR::openMSfile(file)
    temp <- mzR::header(mzRobj)
    #positions to catch
    pos2Catch <- temp$msLevel == msLevel
    if(!missing(acquisitionNum)){
      pos2Catch <- pos2Catch & temp$acquisitionNum %in% acquisitionNum
    }
    #not an scan to read, jump to next file
    if(!any(pos2Catch)) next

    #load metadata & filename. Append filename to spectra metadata
    temp <- cbind(file = basename(file), temp[pos2Catch,])

    # load spectra matrix
    temp_spctra <- lapply(which(pos2Catch), function(x) {
      m <- t(mzR::peaks(mzRobj,x))
      rownames(m) <- c("mass-charge","intensity")
      return(m)
    })

    #merge info
    spctra <- c(temp_spctra, spctra)
    mtdata <- rbind(temp,mtdata)
    mzR::close(mzRobj)
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
      set.seed(1)
        spectra2subset <- sample(mtdata$idSpectra, nsamples)
        mtdata <- mtdata[mtdata$idSpectra %in% spectra2subset, ]
        spctra <- spctra[names(spctra) %in% spectra2subset]
      }
  }
  result <- list(Metadata = mtdata, Spectra = spctra)
  return(result)
}
