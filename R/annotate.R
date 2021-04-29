#' Annotate spectra against a MS2ID library
#'
#' \code{annotate} returns, for every query spectrum, a list of compound
#' candidates from a reference library (\linkS4class{MS2ID} object). The
#' criteria relies on compare every query spectrum with reference spectra using
#' distance metrics.
#'
#' @param QRYdir character(1). Path of the folder containing the mzML files
#'   (query spectra) to annotate.
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
#' @param massError TODO
#' @param cmnPeaks -Reference spectra filter- Integer limiting reference spectra
#'   to whose with at least that number of peaks in common with the query
#'   spectrum.
#' @param cmnTopPeaks -Reference spectra filter- Integer limiting the annotation
#'   to reference spectra with at least one common peak -with the query
#'   spectrum- among their top n most intense peaks.
#' @param cmnPolarity -Reference spectra filter- Boolean, a TRUE value limits
#'   the reference spectra to those with the same polarity as the query
#'   spectrum.
#' @param db -Reference spectra filter- Character filtering the reference
#'   spectra by its original database (e.g. HMDB).
#' @param nature -Reference spectra filter- Character filtering the reference
#'   spectra by the spectra nature. Values are restricted to 'experimental',
#'   'insilico' or 'all'.
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
#' result <- annotate(QRYdir = q, MS2ID = ms2idObject, nsamples=10,
#'          metrics = c("fidelity", "cosine", "topsoe"),
#'          metricsThresh = c(0.6, 0.8, 0.6),
#'          metricFUN = fooFunction, metricFUNThresh = 1.8)
#' }

annotate <- function(QRYdir, MS2ID, metrics="cosine", metricsThresh= 0.8,
                     metricFUN, metricFUNThresh, noiseThresh=0.01, massError=20,
                     cmnPrecMass=FALSE, cmnNeutralMass=TRUE, cmnPeaks=2,
                     cmnTopPeaks=5, cmnPolarity= TRUE, db="all", nature="all",
                     nsamples, consens=T, consCos=0.8, consComm=2/3, consME=20,
                     ...){
    argmnts <- c(as.list(environment()), list(...))
    #check argument types
    reqClasses <- c(QRYdir="character", MS2ID="MS2ID", metricsThresh="numeric",
                    metricFUNThresh="numeric", nsamples="integer",
                    noiseThresh="numeric", cmnPeaks="integer",
                    cmnTopPeaks="integer", db="character", nature="character",
                    cmnPolarity= "logical", cmnPrecMass="logical",
                    cmnNeutralMass="logical", massError="numeric")

    reqClasses <- reqClasses[names(reqClasses) %in% names(argmnts)]
    .checkTypes(argmnts[match(names(reqClasses), names(argmnts))], reqClasses)

    #type of metric (incrm. or decremental)
    decrMetrics <- c("topsoe", "squared_chord")
    decrMet <- metrics %in% decrMetrics

    #check values
    if(!all(metrics %in% c("cosine", "topsoe", "fidelity", "squared_chord")))
      stop(paste("'metrics' is expected to be one of the following characters:",
                 "'cosine', 'topsoe', 'fidelity', 'squared_chord'"))
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

    if (missing(QRYdir))
        stop("'QRYdir' is a mandatory argument")
    if(noiseThresh < 0 | noiseThresh > 1)
        stop("'noiseThresh' is expected to be a number between 0 and 1")
    if (cmnPeaks < 1)
        stop("'cmnPeaks' is expected to be a natural number")
    if (!cmnTopPeaks %in% 1:6 )
        stop("'cmnTopPeaks' is expected to be a natural number between 1 and 6")
    if (!nature %in% c("experimental", "insilico", "all"))
        stop("'nature' is expected to be 'experimental', 'insilico' or 'all'")
    if(massError < 0)
        stop("'massError' is expected to be a positive number")

    workVar <- mget(names(formals()), sys.frame(sys.nframe()))
    workVar$MS2ID <- dirname(MS2ID@dbcon@dbname)
    workVar$annotationTime <- Sys.time()

    message("Loading query spectra ...")
    QRY <- .loadSpectra(dirPath=QRYdir, nsamples=nsamples, ...)

    #check queryif argument "cmnPolarity" can be applied
    if(cmnPolarity){
        if(!all(QRY$Metadata$polarity %in% c(0,1)))
            stop("cmnPolarity=TRUE can not be applied because some query spectra
                 have unknown polarity (neither 1 (positive) nor 0 (negative))")
    }
    #remove invalid spectra
    QRY <- .validateSpectra(QRY)

    # Clean spectra. Remove fragments with intensity < 1% base peak
    # (considering noiseThresh=0.01)
    QRY$Spectra <- lapply(QRY$Spectra, function(x) {
        x[,x["intensity",] > noiseThresh * max(x["intensity", ]), drop = F]
    })

    # bin spectra
    #TODO: Make this cte a variable (dependent of the mz index resolution) and
    # check if its consistent with MS2ID bin
    dec2binFrag=2
    message("Binning query spectra ...")
    QRY$Spectra <- .binSpectra(spectraList=QRY$Spectra, dec2binFrag)

    if(consens){
      message("Obtaining consensus spectra ...")
      #cluster spectra to consens them
      QRY <- .cluster2Consens(QRY, consCos, massErrorPrec = consME)
      #consens spectra
      QRY <- .consens(QRY, consME, consComm)
      #keep only consensus spectra (to annotate them)
      isConsSpectra <- vapply(QRY$Metadata$primalSpctra, function(primals){
        !is.null(primals) & !identical(primals,0)
      }, FUN.VALUE = T)
      if(sum(isConsSpectra)==0) stop("Not even one consensus spectrum obtained")
      QRY$Metadata <- QRY$Metadata[isConsSpectra,]
      QRY$Spectra <- QRY$Spectra[names(QRY$Spectra) %in% QRY$Metadata$idSpectra]
      QRY$Spectra <- .binSpectra(spectraList=QRY$Spectra, dec2binFrag)
    }

    # SQL sentence according global restrictions (db, nature, ...
    SQLwhereGen <- vector()
    if(db!="all")
        SQLwhereGen <- .appendSQLwhere("REFID_db", db,
                                       whereVector=SQLwhereGen)
    if(nature!="all")
        SQLwhereGen <- .appendSQLwhere("REFnature", nature,
                                       whereVector = SQLwhereGen)
    message("Solving distance metrics between query and reference spectra ...")

    distances <- pbapply::pblapply(seq_along(QRY$Spectra), function(idQspctr){
        posMetadata <- which(QRY$Metadata$idSpectra ==
                                 names(QRY$Spectra[idQspctr]))

        massQ <- QRY$Spectra[[idQspctr]]["mass-charge",]
        intQ <- QRY$Spectra[[idQspctr]]["intensity",]

        idRef <- .queryMzIndex(mzVector = massQ, ms2idObj = MS2ID,
                               cmnPeaks = cmnPeaks, cmnTopPeaks = cmnTopPeaks)

        SQLwhereIndv <- vector()
        #apply polarity filter
        if(cmnPolarity){
            SQLwhereIndv <- .appendSQLwhere("REFpolarity",
                                            QRY$Metadata$polarity[posMetadata],
                                            whereVector=SQLwhereIndv)
        }
        #apply precursor mass filter
        if(cmnPrecMass){
            minPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 - massError/10^6)
            maxPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 + massError/10^6)
            SQLwhereIndv <- .appendSQLwhere("REFprecursor_mz",
                                            c(minPrecMass, maxPrecMass),
                                            mode="BETWEEN",
                                            whereVector=SQLwhereIndv)
        }

        if(cmnNeutralMass){
          QRYMmi <- .propQMmi(QRY$Metadata$precursorMZ[posMetadata],
                              QRY$Metadata$polarity[posMetadata])
          QRYMmi_min <- QRYMmi * (1 - massError/10^6)
          QRYMmi_max <- QRYMmi * (1 + massError/10^6)

          subSQLwhereMmi <- .appendSQLwhere("REFMmi",
                                            c(min(QRYMmi_min), max(QRYMmi_max)),
                                            mode="BETWEEN")
          idRefMmi <- .getSQLrecords(MS2ID, select="ID_metabolite, REFMmi",
                                     from="metaCompound",where=subSQLwhereMmi)
          idRefMeta <- idRefMmi$ID_metabolite[vapply(idRefMmi$REFMmi, function(ix)
            any(ix > QRYMmi_min & ix < QRYMmi_max), FUN.VALUE = T)]
          #idRefMeta <- idRefMmi$ID_metabolite[idRefMmi$REFMmi %in% REFMmi_r]
          subSQLwhereMmi <- .appendSQLwhere("ID_metabolite", idRefMeta,
                                            mode="IN")
          #subSQLwhereMmi <- .appendSQLwhere("ID_metabolite", "234263", mode="IN")
          idRefMmi <- .getSQLrecords(MS2ID, select="ID_spectra",
                                     from="crossRef_SpectrComp",
                                     where=subSQLwhereMmi)
          idRef <- intersect(idRef, unlist(idRefMmi))
        }
        SQLwhereIndv <- .appendSQLwhere("ID_spectra", idRef, mode="IN",
                                        whereVector=SQLwhereIndv)
        idRef <- .getSQLrecords(MS2ID,"ID_spectra", "metaSpectrum",
                                c(SQLwhereGen, SQLwhereIndv))

        #return if query spectrum has no targeted db spectra
        if(nrow(idRef) == 0) return(NA)

        #get spectra from big memory
        refSpectra <- .bufferSpectra(MS2ID, idRef$ID_spectra)

        distance <- lapply(seq_along(refSpectra$ptr$id), function(x) {
            a <- .getSpectrum(refSpectra, x)
            row1<- unique(c(a["mass-charge",], massQ))
            row2 <- a["intensity",][match(row1, a["mass-charge",])]
            row1 <- intQ[match(row1, massQ)]
            row1[is.na(row1)] <- 0
            row2[is.na(row2)] <- 0
            #normalize intensities and add 1e-12 (2 avoid problems with log(0))
            rowdf <- rbind(row1/sum(row1), row2/sum(row2)) + 1e-12
            rsltM <- vapply(metrics, function(iM){
              suppressMessages(philentropy::distance(rowdf, method = iM))
            }, FUN.VALUE = 3.2)
            if(metFun) rsltM <- c(rsltM, metricFunc =
                                    metricFUN(rbind(massQ, intQ), a))
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
          distance <- distance[hits, , drop=F]
          distance$idREFspect <- refSpectra$ptr$id[hits]
          return(distance)
        } else {return(NA)}

    })

    names(distances) <- names(QRY$Spectra)
    #remove query spectra with no hits
    distances <- distances[!is.na(distances)]

    if(length(distances) == 0){
        stop("No query spectra with hits found.")
    }
    message("Processing results obtained ...")
    rslt <- .processAnnotation(dist = distances, ms2id = MS2ID , qry = QRY,
                               mError = massError, cmnNtMass = cmnNeutralMass,
                               workVar = workVar)
    return(rslt)
}
#TODO: SHINY
#TODO: CONSENSUS SPECTRA
