#'
#' Annotate MS/MS spectra against a reference library
#'
#' \code{annotate} is the MS2ID function that annotates MS/MS query spectra.
#' Every query spectrum is compared with a reference library
#' (\linkS4class{MS2ID} object) using as criteria different distance metrics;
#' the function returns query and reference spectra (and their compounds) that
#' beat a determined threshold.
#'
#' @param QRYdata Query spectra. \code{\link[Spectra]{Spectra}}(1) object (see
#'   \href{https://bioconductor.org/packages/3.14/Spectra}{Spectra} package) or
#'   character(1) with the directory name containing the mzML files.
#' @param MS2ID \linkS4class{MS2ID} object with the in-house database to use in
#'   the annotation.
#' @param noiseThresh numeric(1) defining the threshold used in the noise
#'   filtering of the query spectra. It is expressed as \% intensity relative to
#'   base peak. e.g. noiseThresh=0.01 eliminates peaks with an intensity of less
#'   than 1\% of the base peak.
#' @param metrics character(n) that defines which metrics use to compare query
#'   and reference spectra (\code{annotate} function trims the result according
#'   the \code{metricsThresh} parameter). Values are restricted to 'cosine',
#'   'topsoe', 'fidelity' and 'squared_chord'. See
#'   \code{\link[philentropy]{distance}} function in philentropy package for
#'   more information.
#' @param metricsThresh numeric(n) defining the n threshold values of the n
#'   metrics. A reference spectrum is considered a hit and listed in the return
#'   object when \strong{at least one of the metrics} fulfills its threshold;
#'   note that to fulfill a threshold value has a different meaning depending on
#'   the metric: topsoe and squared_chord metrics return a lower number when the
#'   spectra are more similar so, unlike the rest, a hit will occur when the
#'   returned value is lower than its threshold. Recommended values to start
#'   with are cosine=0.8, topsoe=0.6, fidelity=0.6 and squared_chord=0.8.
#' @param metricFUN function(1) defined by the user to be used as metric. This
#'   function must accept a two-row matrix as a parameter; each row must contain
#'   the intensity of a spectrum, and each column must refer to the same m/z
#'   value (considering the massErrMsn parameter as the mass error). Finally,
#'   the function must return a numeric(1). See example.
#' @param metricFUNThresh numeric(1) with the threshold value of the metric
#'   defined by the \code{metricFUN} parameter. metricFUN / metricFUN are
#'   analogous to metrics / metricsThresh parameters.
#' @param massErrMs1 numeric(1). Mass error to consider in operations with first
#'   spectrometer measures (MS1), e.g. grouping spectra according its precursor
#'   mass (in consensus formation) or evaluating precursor and neutral masses
#'   similarities (in reference spectra prefiltering),
#' @param massErrMsn numeric(1) Mass error to consider in operations with
#'   non-first spectrometer measures (typically MS2). e.g. matching fragments,
#'   for consensus formation or distance similarity measures.
#' @param cmnFrags -Reference spectra filter- vector with two integers (m, n)
#'   limiting the reference spectra so that both query and reference spectra
#'   have at least m peaks in common among their top n most intense peaks.
#' @param cmnPolarity -Reference spectra filter- Boolean(1) that limits the
#'   reference spectra to those with the same polarity than the query spectrum.
#' @param predicted -Reference spectra filter- Boolean(1) filtering the
#'   reference spectra according its experimental nature (in-silico or
#'   predicted). A NULL value does not apply filter.
#' @param cmnPrecMass -Reference spectra filter- Boolean(1) that limits the
#'   reference spectra to those that have the precursor mass of the query
#'   spectrum.
#' @param cmnNeutralMass -Reference spectra filter- Boolean filtering the
#'   reference spectra to those with a neutral mass plausible with the query
#'   precursor (considering all possible adducts).
#' @param nsamples integer(1) defines a subset of x random query spectra to work
#'   with. Useful for speeding up preliminary testing before definitive
#'   annotation, it is not compatible with the consensus formation.
#'
#' @return an \linkS4class{Annot} object with the results of the annotation
#' @export
#' @example man/examples/annotate.R
annotate <- function(QRYdata, QRYmsLevel = 2L, MS2ID,
                     metrics="cosine", metricsThresh= 0.8,
                     metricFUN, metricFUNThresh,
                     massErrMs1 = 5, massErrMsn = 20,
                     noiseThresh = 0.01,  cmnPrecMass = FALSE,
                     cmnNeutralMass = TRUE, cmnFrags = c(2,5),
                     cmnPolarity = TRUE, predicted = NULL,
                     consens=TRUE, consCos=0.8, consComm=2/3,
                     ...){
  argmnts <- c(as.list(environment()), list(...))
  if(is.na(QRYmsLevel)) QRYmsLevel <- NULL
  if(length(cmnFrags)==2 & is.numeric(cmnFrags)){
    if(cmnFrags[1] < 1 | cmnFrags[1] != as.integer(cmnFrags[1]))
      stop(glue::glue("First position of 'cmnFrags' argument is expected to \\
                      be a integer > 1"))
    if(!cmnFrags[2] %in% 1:6)
      stop(glue::glue("Second position of 'cmnFrags' argument is expected to \\
                      be a integer between 1 and 6"))
    if(!cmnFrags[2] %in% 1:6)
      stop(glue::glue("In the 'cmnFrags' argument, first integer can not be \\
      greater than the second one. e.g. it is not possible to find spectra \\
      with 4 peaks in common among its top 3 most intense peaks
                      "))
  }else
    stop("'cmnFrags' argument is expected to be a vector of 2 integers")

    #check argument types
  reqClasses <- c(MS2ID="MS2ID", QRYmsLevel = "integer",
                  metricsThresh="numeric", metricFUNThresh="numeric",
                  noiseThresh="numeric", predicted="logical",
                  cmnPolarity= "logical", cmnPrecMass= "logical",
                  cmnNeutralMass="logical",
                  massErrMs1="numeric", massErrMsn="numeric")

    reqClasses <- reqClasses[names(reqClasses) %in% names(argmnts)]
    .checkTypes(argmnts[match(names(reqClasses), names(argmnts))], reqClasses)

    #type of metric (incrm. or decremental)
    decrMet <- metrics %in% DECRMETRIC

    metOK <- c(INCRMETRIC,DECRMETRIC)[c(INCRMETRIC,DECRMETRIC)!='metricFunc']
    #check values
    if(!all(metrics %in% metOK))
      stop(paste("'metrics' is expected to be one of the following:",
                 paste(metOK, collapse = ", ")))
    if(length(metrics) != length(metricsThresh))
      stop("'metricsThresh' must contain a value for every metric included
             in the argument 'metrics'")
    if(!missing(metricFUN)){
      if(!is.function(metricFUN))
        stop("'metricFUN' argument is expected to be a function")
      if(missing(metricFUNThresh))
        stop(paste("A numeric is expected in 'metricFUNThresh' when metricFUN",
                   "argument is used is expected to a numeric "))
      metFun <- TRUE
      metricsThresh <- c(metricsThresh, metricFUNThresh)
      decrMet <- c(decrMet, F)
      } else metFun <- FALSE

    if (missing(QRYdata))
        stop("'QRYdata' is a mandatory argument")
    if(noiseThresh < 0 | noiseThresh > 1)
        stop("'noiseThresh' is expected to be a number between 0 and 1")
    if(massErrMs1 < 0)
      stop("'massErrMs1' is expected to be a positive number")
    if(massErrMsn < 0)
      stop("'massErrMsn' is expected to be a positive number")

    workVar <- mget(names(formals()), sys.frame(sys.nframe()))
    if(!is.character(workVar$QRYdata)) workVar$QRYdata <- "spectra object"
    workVar$MS2ID <- dirname(MS2ID@dbcon@dbname)
    workVar$annotationTime <- Sys.time()

    message("Loading query spectra ...")
    QRY <- .loadSpectra(data = QRYdata, msLevel = QRYmsLevel, ...)
    #Check arguments' viability considering query metadata availability
    if(consens)
      .checkViability(metadata=QRY$Metadata, argument="consens",
                      necMetadata=c("precursorMZ", "precursorIntensity",
                                    "collisionEnergy","file", "polarity"))
    if(cmnPolarity){
      .checkViability(metadata=QRY$Metadata, argument="cmnPolarity",
                      necMetadata=c("polarity"))
      if(!all(QRY$Metadata$polarity %in% c(0, 1)))
        stop(glue::glue("
            cmnPolarity = TRUE can not be applied because some query spectra \\
            have unknown polarity (neither 1 (positive) nor 0 (negative))
                            "))
    }
    if(cmnPrecMass)
      .checkViability(metadata=QRY$Metadata, argument="cmnPrecMass",
                      necMetadata=c("precursorMZ"))
    if(cmnNeutralMass)
      .checkViability(metadata=QRY$Metadata, argument="cmnNeutralMass",
                      necMetadata=c("precursorMZ", "polarity"))
    #remove invalid spectra
    QRY <- .validateSpectra(QRY)

    # Clean spectra. Remove fragments with intensity < 1% base peak
    # (considering noiseThresh=0.01)
    QRY$Spectra <- lapply(QRY$Spectra, function(x) {
        x[,x["intensity",] > noiseThresh * max(x["intensity", ]), drop = F]
    })

    LFT <- NA

    if(consens){
      message("Obtaining consensus spectra ...")
      #cluster spectra to consens them
      QRY <- .cluster2Consens(s = QRY, consCosThres = consCos,
                              massErrMs1 = massErrMs1, massErrMsn = massErrMsn)

      if(all(QRY$Metadata$rol != 4L)) {
        message("No consensus spectrum was obtained")
      }else{
        #consens the spectra
        QRY <- .consens(QRY,  massError = massErrMsn, consComm)
      }
      #LFT: leftovers, spectra not to annotate only to keep query spectra temporaly just for traceability of consensus formation
      rol2Annotate <- c(1, 2, 4)
      LFT <- QRY
      LFT$Metadata <- QRY$Metadata[!QRY$Metadata$rol %in% rol2Annotate, ]
      LFT$Spectra <- QRY$Spectra[names(QRY$Spectra) %in%
                                   LFT$Metadata$idSpectra]
      QRY$Metadata <- QRY$Metadata[QRY$Metadata$rol %in% rol2Annotate, ]
      QRY$Spectra <- QRY$Spectra[names(QRY$Spectra) %in% QRY$Metadata$idSpectra]
    }else{
      QRY$Metadata$rol <- 1L
    }

    # SQL sentence according global restrictions (predicted)
    SQLwhereGen <- vector()
    if(!is.null(predicted)){
      predicted <- ifelse(predicted, 1, 0)
      SQLwhereGen <- .appendSQLwhere("predicted", predicted,
                                     whereVector = SQLwhereGen)
    }
    message("Solving distance metrics between query and reference spectra ...")

    distances <- pbapply::pblapply(seq_along(QRY$Spectra), function(idQspctr){
      posMetadata <- which(QRY$Metadata$idSpectra ==
                             names(QRY$Spectra[idQspctr]))
      Qspct <- QRY$Spectra[[idQspctr]]
      idRef <- .queryMzIndex(QRYspct = Qspct, ms2idObj = MS2ID,
                             cmnFrags = cmnFrags)
      Qspct <- rbind(Qspct,
                     error = Qspct["mass-charge", ] * massErrMsn/1e6,
                     intSpectr2 = 0)
      SQLwhereIndv <- vector()
        #apply polarity filter
        if(cmnPolarity){
            SQLwhereIndv <- .appendSQLwhere("polarity",
                                            QRY$Metadata$polarity[posMetadata],
                                            whereVector=SQLwhereIndv)
        }
        #apply precursor mass filter
        if(cmnPrecMass){
            minPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 - massErrMs1/1e6)
            maxPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 + massErrMs1/1e6)
            SQLwhereIndv <- .appendSQLwhere("precursorMz",
                                            c(minPrecMass, maxPrecMass),
                                            mode="BETWEEN",
                                            whereVector=SQLwhereIndv)
        }

        if(cmnNeutralMass){
          QRYMmi <- .propQMmi(QRY$Metadata$precursorMZ[posMetadata],
                              QRY$Metadata$polarity[posMetadata])
          QRYMmi_min <- QRYMmi * (1 - massErrMsn/1e6)
          QRYMmi_max <- QRYMmi * (1 + massErrMsn/1e6)
          subSQLwhereMmi <- .appendSQLwhere("exactmass",
                                            c(min(QRYMmi_min), max(QRYMmi_max)),
                                            mode="BETWEEN")
          idRefMmi <- .getSQLrecords(MS2ID,
                                     select="ID_compound, exactmass",
                                     from="metaCompound",
                                     where = subSQLwhereMmi)
          idRefMeta <- idRefMmi$ID_compound[vapply(idRefMmi$exactmass,
                                                   function(ix)
            any(ix > QRYMmi_min & ix < QRYMmi_max), FUN.VALUE = T)]
          subSQLwhereMmi <- .appendSQLwhere("ID_compound", idRefMeta,
                                            mode="IN")
          idRefMmi <- .getSQLrecords(MS2ID, select="ID_spectra",
                                     from = "crossRef_SpectrComp",
                                     where = subSQLwhereMmi)
          idRef <- intersect(idRef, unlist(idRefMmi))
        }
        SQLwhereIndv <- .appendSQLwhere("ID_spectra", idRef, mode="IN",
                                        whereVector=SQLwhereIndv)
        idRef <- .getSQLrecords(MS2ID, "ID_spectra", "metaSpectrum",
                                c(SQLwhereGen, SQLwhereIndv))

        #return if query spectrum has no targeted db spectra
        if(nrow(idRef) == 0) return(NA)

        #get spectra from big memory
        refSpectra <- .bufferSpectra(MS2ID, idRef$ID_spectra)
        distance <- lapply(seq_along(refSpectra$ptr$id), function(x) {
          struct <- .matchFrag(Qspct, .getSpectrum(refSpectra, x))
          #normalize intensities and add 1e-12 (2 avoid problems with log(0))
          rowdf <- rbind(struct[1,]/sum(struct[1,]),
                         struct[2,]/sum(struct[2,])) + 1e-12
          rsltM <- vapply(metrics, function(iM){
            suppressMessages(philentropy::distance(rowdf, method = iM))
          }, FUN.VALUE = 3.2)
          if(metFun){
            rsltM <- c(rsltM, metricFunc = metricFUN(rowdf))
          }
          return(rsltM)
        })
        distance <- data.frame(do.call(rbind, distance))
        hits <- lapply(seq_len(ncol(distance)), function(im){
          if(decrMet[im]) distance[, im] <= metricsThresh[im]
          else  distance[, im] >= metricsThresh[im]
        })
        hits <- do.call(cbind, hits)
        hits <- apply(hits, MARGIN = 1, any)
        if(any(hits)){
          distance <- distance[hits, , drop = F]
          distance$idREFspect <- refSpectra$ptr$id[hits]
          return(distance)
        } else {return(NA)}
    })

    names(distances) <- names(QRY$Spectra)
    #remove query spectra with no hits
    distances <- distances[!is.na(distances)]
    if("rawSpectra" %in% names(QRY)){
      QRY$Spectra <- QRY$rawSpectra
      QRY$rawSpectra <- NULL
    }

    if(length(distances) == 0){
        message("No query spectrum has obtained a hit")
        return()
    }
    message("Processing results obtained ...")
    rslt <- .processAnnotation(dist = distances, ms2id = MS2ID ,
                               qry = QRY, lft = LFT,
                               massErrMs1 = massErrMs1, massErrMsn = massErrMsn,
                               cmnNtMass = cmnNeutralMass, workVar = workVar)
    return(rslt)
}

.checkViability <- function(metadata, argument, necMetadata){
  matchnecCols <- necMetadata %in% names(metadata)
  if(!all(matchnecCols) | anyNA(metadata[, necMetadata[matchnecCols]]))
    stop(glue::glue("
      '{argument}' can NOT be applied: query spectra have incomplete some of \\
      the following required metadata: \\
      {glue::glue_collapse(necMetadata, ', ', last = ' or ')}
                      "))
}
