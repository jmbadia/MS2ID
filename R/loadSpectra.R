#' Load query spectra
#'
#' @param data `character` defining either the directory containing the mzML files, or the files themselves.
#' @param msLevel `integer` or `numeric`, but a natural number. It subsets the
#'  spectra to be loaded to those with such msLevel. By default, msLevel=2L
#' @param nsamples integer(1) defines a subset of x random query spectra to work
#'   with. Useful for speeding up preliminary testing before definitive
#'   annotation, it is not compatible with the consensus formation.
#'
#' @param acquisitionNum vector of`integer` or `numeric`. The user can limit spectra based on their acquisitionNum
#'
#' @return a list with 2 items, 'Metadata' and 'Spectra'. The former is a data frame with spectrum metadata. The latter is a list with two items, a list of spectra (under matrix form) and 'idSpectra' (vector of its spectra id.) Both 'Metadata' and 'Spectra' are linked using the 'idSpectra' variable.
#' @noRd

.loadSpectra <- function(data = NULL, nsamples = NULL, ...){
  argmnts <- c(as.list(environment()), list(...))
  if(is.null(data))
    stop("'data' is a mandatory argument")
  #check types
  reqClasses <- c(nsamples = "integer")
  .checkTypes(argmnts, reqClasses)
  if (!is.null(nsamples))
    if (nsamples < 1)
      stop("'nsamples' is expected to be a natural number > 0")
  if(is.character(data)){#if data are path to files
    if(identical(file_test("-d", data), T)){# one char and is a dir
      #load all valid files
      data <- dir(path = data, pattern="*\\.mzml$|*\\.msp$|*\\.mgf$",
                       ignore.case = TRUE, full.names = TRUE)
    }
    ismzml <- grepl("\\.mzml$", data, ignore.case = TRUE)
    ismsp <- grepl("\\.msp$", data, ignore.case = TRUE)
    ismgf <- grepl("\\.mgf$", data, ignore.case = TRUE)
    formats <- c("mzml", "msp", "mgf")
    formatPresence <- c(any(ismzml), any(ismsp), any(ismgf))
    if(sum(formatPresence) > 1){
      formats <- c("mzml", "msp", "mgf")[formatPresence]
      stop(glue::glue("Detected \\
        {glue::glue_collapse(formats, ', ', last = ' and ')} file \\
        formats. Only one type of file format is simultaneously allowed."))
    }else if(sum(formatPresence) < 1)
      stop(glue::glue("No valid file formats \\
        ({glue::glue_collapse(formats, ', ', last = ' or ')}) detected."))
    validFile <- ismzml | ismsp | ismgf
    if(!all(validFile))
      warning(paste("The following non valid file/s will be ignored:",
                    paste(data[!validFile], collapse=", ")))
    if(any(ismzml))
      result <- .loadmzML(mzmlData = data[ismzml], ...)
    if(any(ismsp)){
      bckend <- backendInitialize(MsBackendMsp::MsBackendMsp(), data[ismsp])
      spectra <- Spectra::Spectra(bckend)
      result <- .convertSpectra(spectra = spectra, ...)
    }
    if(any(ismgf)){
      bckend <- backendInitialize(MsBackendMgf::MsBackendMgf(), data[ismgf])
      spectra <- Spectra::Spectra(bckend)
      result <- .convertSpectra(spectra = spectra, ...)
    }
  }else if(is(data, 'Spectra')){
    result <- .convertSpectra(spectra = data, ...)
  } else {
    stop(glue::glue("'data' is expected to be an spectra object or a //
      character describing a directory's path with mzML files"))
  }
# jmbadia comprovar loadmzml no te redundancies respecte la funcio mare. I ficar llibrarries
#subset to a number 'nsamples' of random samples.
# In this part of the code in order to maintain id traceability
if(!is.null(nsamples)){
  if(nsamples > nrow(result$Metadata))
      message(glue::glue("
      'nsamples' is larger than the number of spectra available. All //
      spectra will be loaded"))
  else{
      spectra2subset <- sample(result$Metadata$idSpectra, nsamples)
      result$Metadata <- result$Metadata[
        result$Metadata$idSpectra %in% spectra2subset, ]
      result$Spectra <- result$Spectra[
        names(result$Spectra) %in% spectra2subset]
    }
}
return(result)
}

#' Convert spectra (Spectra package) to MS2ID query format
#'
#' @noRd
.convertSpectra <- function(spectra=NULL, msLevel=NULL, acquisitionNum=NULL){
  if(is.null(spectra))
    stop("'spectra' is a mandatory argument")
  #check types
  reqClasses <- c(spectra="Spectra", msLevel="integer",
                  acquisitionNum = "integer")
  argmnts <- c(as.list(environment()))
  .checkTypes(argmnts, reqClasses)

  if(!is.null(msLevel))
    if (msLevel < 1)
      stop("'msLevel' is expected to be a natural number > 0")
  if (!is.null(acquisitionNum))
    if (any(acquisitionNum < 1))
      stop("'acquisitionNum' are expected to be natural numbers > 0")
  if(!is.null(msLevel))
    spectra <- Spectra::filterMsLevel(spectra, msLevel = msLevel)
  else
    msLevel <- "all"
  if(!is.null(acquisitionNum))
    spectra <- Spectra::filterAcquisitionNum(spectra, acquisitionNum)
  else
    acquisitionNum <- "all"
  if(length(spectra) == 0)
    stop(glue::glue("
    No query spectra match the arguments used in the \\
    annotate() function (msLevel={msLevel}, acquisitionNum={acquisitionNum})
                    "))
  #rename only if it exists
  possibleCols <-
    list(
      precursorScanNum = "precScanNum",
      isolationWindowTargetMZ = "isolationWindowTargetMz"
    )
  possibleCols <- possibleCols[possibleCols %in% intersect(
    unlist(possibleCols),
    Spectra::spectraVariables(spectra)
  )]
  mtdata <- spectra %>% Spectra::spectraData() %>%
    as.data.frame() %>%
    rename(retentionTime = rtime, file = dataOrigin,
           precursorMZ = precursorMz, !!!possibleCols) %>%
    select(!dataStorage)
  mtdata$file <- basename(mtdata$file)
  spctra <- lapply(Spectra::peaksData(spectra), function(idxsp) {
    colnames(idxsp) <- c('mass-charge', 'intensity')
    t(idxsp)
  })

  # Apply same ID to metadata spectra and metadata matrix
  mtdata <- cbind(idSpectra = seq_len(nrow(mtdata)), mtdata)
  names(spctra) <- seq_len(nrow(mtdata))

  return(list(Metadata = mtdata, Spectra = spctra))
}

#' Load query spectra from mzML files
#'
#' @noRd
.loadmzML <- function(mzmlData=NULL, msLevel=NULL, acquisitionNum=NULL){
  if(is.null(mzmlData))
    stop("'mzmlData' is a mandatory argument")
  reqClasses <- c(mzmlData="character",
                  msLevel="integer", acquisitionNum = "integer")
  .checkTypes(c(as.list(environment())), reqClasses)
  if (msLevel < 1)
    stop("'msLevel' is expected to be a natural number > 0")
  if (!is.null(acquisitionNum))
    if (any(acquisitionNum < 1))
      stop("'acquisitionNum' are expected to be natural numbers > 0")
  ismzml <- grepl("\\.mzml$", mzmlData, ignore.case = TRUE)
  if(!all(ismzml))
    warning(paste("The following non mzML file/s will be ignored:",
                  paste(mzmlData[!ismzml], collapse=", ")))
  mzml_files <- mzmlData[ismzml]
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
    pos2Catch <- rep(TRUE, nrow(temp))
    if(!is.null(msLevel))
      pos2Catch <- pos2Catch & temp$msLevel == msLevel
    if(!is.null(acquisitionNum))
      pos2Catch <- pos2Catch & temp$acquisitionNum %in% acquisitionNum
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
    mtdata <- rbind(temp, mtdata)
    mzR::close(mzRobj)
  }

  # Apply same ID to metadata spectra and metadata matrix
  mtdata <- cbind(idSpectra = seq_len(nrow(mtdata)), mtdata)
  names(spctra) <- seq_len(nrow(mtdata))
  return(list(Metadata = mtdata, Spectra = spctra))
}
