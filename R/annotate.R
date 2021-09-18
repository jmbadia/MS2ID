#'
#' Annotate spectra against a MS2ID library
#'
#' \code{annotate} returns, for every query spectrum, a list of compound
#' candidates from a reference library (\linkS4class{MS2ID} object). The
#' criteria relies on compare every query spectrum with reference spectra using
#' distance metrics.
#'
#' @param QRYdata character(n) defining either the directory containing the mzML files, or the files themselves.
#' @param MS2ID \linkS4class{MS2ID} object with the in-house database to use in
#'   the annotation.
#' @param noiseThresh A numeric defining the threshold used in the noise
#'   filtering of the query spectra, considered as \% intensity relative to base
#'   peak. e.g. noiseThresh=0.01 eliminates peaks with an intensity of less than
#'   1\% of the base peak.
#' @param metrics character(n) defining the n distance metrics to measure
#'   simultaneously between query and reference spectra. Values are restricted
#'   to 'cosine', 'topsoe', 'fidelity' and 'squared_chord'.
#' @param metricsThresh numeric(n) defining the n threshold values of the n
#'   metrics. A reference spectrum is considered a hit when at least one of the
#'   metrics measured fulfills its threshold. Recommended values to start with
#'   are cosine=0.8, topsoe=0.6, fidelity=0.6 and squared_chord=0.8.
#' @param metricFUN function(1) user-made defining a new metric for annotation.
#'   The function must have the two spectra to be compared as arguments. Each of
#'   these spectra must be a matrix with its mass-charge and intensity in rows
#'   respectively. Finally, the function must return a numeric(1).
#' @param metricFUNThresh numeric(1) threshold value of the metric defined by
#'   the 'metricFUN' function.
#' @param massErrMs1 numeric(1). Mass error to consider in operations with first spectrometer measures, e.g.  in consensus formation, grouping query spectra according its precursor mass or, in reference spectra prefiltering, evaluate its precursor or neutral masses
#' @param massErrMsn numeric(1) Mass error to consider in operations with non first spectrometer measures (typycally MS2). e.g. grouping fragments for consensus formation or distance similarity measures (tipically cosine)
#'
#' @param cmnFrags -Reference spectra filter- vector with 2 integers (m, n) limiting reference spectra. Reference spectra and query spectra must have at least m peaks in common among their top n most intense peaks
#' @param cmnPolarity -Reference spectra filter- Boolean, a TRUE value limits
#'   the reference spectra to those with the same polarity as the query
#'   spectrum.
#' @param predicted -Reference spectra filter- Character filtering the reference
#'   spectra by the spectra nature. Default is no filtering
#' @param cmnPrecMass -Reference spectra filter- Boolean, a TRUE value limits
#'   the reference spectra to those that have the same precursor mass as the
#'   query spectrum.
#' @param cmnNeutralMass -Reference spectra filter- Boolean, a TRUE value limits
#'   the reference spectra to those that have a neutral mass that matches some
#'   of the plausible query neutral masses (considering the precursor query mass
#'   and all the possible adducts, TODO: see link).
#' @param nsamples integer(1) defines a subset of x random query spectra to work
#'   with. Useful for speeding up preliminary testing before definitive
#'   annotation.
#'
#' @return an \linkS4class{Annot} object with the results of the
#'   annotation
#' @export
#' @examples
#' \dontrun{
#' fooFunction <- function(spectr1, spectr2){
#'   mz1 <- spectr1[1, ]
#'   int1 <- spectr1[2, ]
#'   mz2 <- spectr2[1, ]
#'   int2 <- spectr2[2, ]
#'   row1<- unique(c(mz2, mz1))
#'   row2 <- int2[match(row1, mz2)]
#'   row1 <- int1[match(row1, mz1)]
#'   row1[is.na(row1)] <- 0
#'   row2[is.na(row2)] <- 0
#'   rowdf <- rbind(row1, row2)
#'   fooCos <- suppressMessages(philentropy::distance(rowdf, method = "cosine"))
#'   return(fooCos+1)
#' }
#' result <- annotate(QRYdata = q, MS2ID = ms2idObject, nsamples=10,
#'          metrics = c("fidelity", "cosine", "topsoe"),
#'          metricsThresh = c(0.6, 0.8, 0.6),
#'          metricFUN = fooFunction, metricFUNThresh = 1.8)
#' }
annotate <- function(QRYdata, QRYmsLevel = 2L, MS2ID,
                     metrics="cosine", metricsThresh= 0.8,
                     metricFUN, metricFUNThresh,
                     massErrMs1 = 30, massErrMsn = 20,
                     noiseThresh = 0.01,  cmnPrecMass = FALSE,
                     cmnNeutralMass = TRUE, cmnFrags = c(2,5),
                     cmnPolarity = TRUE, predicted=NULL,
                     consens=T, consCos=0.8, consComm=2/3,
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
            rsltM <- c(
              rsltM,
              metricFunc = metricFUN(Qspct[c("mass-charge","intensity"), ], a))
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
