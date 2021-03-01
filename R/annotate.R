#' Annotate spectra against an in-house database
#'
#' \code{annotate} returns, for every query spectrum, a list of compound
#'  candidates from a database by using different metrics to compare their
#'  spectra.
#'
#' @param QRYdir character(1). Path of the folder containing the mzML files
#'  (query spectra) to annotate.
#' @param MS2ID MS2ID object with the in-house database to use in the
#'  annotation.
#' @param noiseThresh A numeric defining the threshold used in the noise
#'  filtering of the query spectra, considered as \% intensity relative to
#'  base peak. e.g. noiseThresh=0.01 eliminates peaks with an intensity of
#'  less than 1\% of the base peak.
#' @param cosSimThresh Numeric defining the minimum cosine similarity a
#'  pair of spectra (query and reference) must have to be returned by the
#'  function. (TODO: see link with the definition of cossim used)
#' @param massError TODO
#' @param cmnPeaks -Reference spectra filter- Integer limiting reference
#' spectra to whose with at least that number of peaks in common with the
#' query spectrum.
#' @param cmnTopPeaks -Reference spectra filter- Integer limiting the
#'  annotation to reference spectra with at least one common peak -with the
#'  query spectrum- among their top n most intense peaks.
#' @param cmnPolarity -Reference spectra filter- Boolean, a TRUE value limits
#'  the reference spectra to those with the same polarity as the query spectrum.
#' @param db -Reference spectra filter- Character filtering the reference
#' spectra by its original database (e.g. HMDB).
#' @param nature -Reference spectra filter- Character filtering the reference
#'  spectra by the spectra nature. Values are restricted to 'experimental',
#'  'insilico' or 'all'.
#' @param cmnPrecMass -Reference spectra filter- Boolean, a TRUE value limits
#'  the reference spectra to those that have the same precursor mass as the
#'  query spectrum.
#' @param cmnNeutralMass -Reference spectra filter- Boolean, a TRUE value
#'  limits the reference spectra to those that have a neutral mass that matches
#'  some of the plausible query neutral masses (considering the precursor
#'  query mass and all the possible adducts, TODO: see link).
#' @param ... Arguments to be passed to internal functions
#'
#' @return
#' @export
#'
#' @examples
#' m <- "/home/jmbadia/Insync/OneDrive/Projectes_R_actius/MS2ID/invisible2Git/resources/MS2ID_B2R_20201113_083214"
#' ms2idObject <- MS2ID(m)
#' q <- "../ID11A19_MS2identification/invisible2Git/MS2ID_Identifications/20201117_MBermudoPOS/rawData"
#' annotate(QRYdir=q, MS2ID=ms2idObject, cmnNeutralMass=FALSE, nsamples=10)

annotate <- function(QRYdir, MS2ID, noiseThresh=0.01, cosSimThresh=0.8,
                     massError=20, cmnPrecMass=FALSE, cmnNeutralMass=TRUE,
                     cmnPeaks=2, cmnTopPeaks=5, cmnPolarity= TRUE, db="all",
                     nature="all", nsamples, ...){

    #check argument types
    reqClasses <- c(QRYdir="character", MS2ID="MS2ID",
                    noiseThresh="numeric", cmnPeaks="integer",
                    cmnTopPeaks="integer", db="character", nature="character",
                    cmnPolarity= "logical", cmnPrecMass="logical",
                    cmnNeutralMass="logical", cosSimThresh="numeric",
                    massError="numeric")
    .checkTypes(as.list(match.call(expand.dots=FALSE))[-1], reqClasses)

    #check values
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
    if(cosSimThresh < 0 | cosSimThresh > 1)
        stop("'cosSimThresh' is expected to be a number between 0 and 1")
    if(massError < 0)
        stop("'massError' is expected to be a positive number")

    annotTime = Sys.time()

    workVar <- mget(names(formals()), sys.frame(sys.nframe()))
    workVar$MS2ID <- dirname(MS2ID@dbcon@dbname)
    workVar <- c(annotationTime=as.character(annotTime), workVar)

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
    QRY$Spectra <- .binSpectra(spectra=QRY$Spectra, dec2binFrag)
    # SQL sentence according global restrictions (db, nature, ...
    SQLwhereGen <- vector()
    if(db!="all")
        SQLwhereGen <- .appendSQLwhere("s.REFID_db", db,
                                       whereVector=SQLwhereGen)
    if(nature!="all")
        SQLwhereGen <- .appendSQLwhere("s.REFnature", nature,
                                       whereVector=SQLwhereGen)
    message("Calculating distance metrics between query and
            reference spectra ...")

    distances <- pbapply::pblapply(seq_along(QRY$Spectra), function(idQspctr){
        posMetadata <- which(QRY$Metadata$idSpectra ==
                                 names(QRY$Spectra[idQspctr]))

        massQ <- QRY$Spectra[[idQspctr]]["mass-charge",]
        intQ <- QRY$Spectra[[idQspctr]]["intensity",]

        idRef <- .queryMzIndex(mzVector=massQ, ms2idObj=MS2ID, cmnPeaks=cmnPeaks,
                               cmnTopPeaks=cmnTopPeaks)

        SQLwhereIndv <- vector()
        #apply polarity filter
        if(cmnPolarity){
            SQLwhereIndv <- .appendSQLwhere("s.REFpolarity",
                                            QRY$Metadata$polarity[posMetadata],
                                            whereVector=SQLwhereIndv)
        }
        #apply precursor mass filter
        if(cmnPrecMass){
            minPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 - massError/10^6)
            maxPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 + massError/10^6)
            SQLwhereIndv <- .appendSQLwhere("s.REFprecursor_mz",
                                            c(minPrecMass, maxPrecMass),
                                            mode="BETWEEN",
                                            whereVector=SQLwhereIndv)
        }

        SQLwhereIndv <- .appendSQLwhere("ID_spectra", idRef, mode="IN",
                                        whereVector=SQLwhereIndv)
        idRef <- .getSQLrecords(MS2ID,"s.ID_spectra", "metaSpectrum s",
                                c(SQLwhereGen, SQLwhereIndv))

        #return if query spectrum has no targeted db spectra
        if(nrow(idRef) == 0) return(NA)

        #get spectra from big memory
        refSpectra <- .bufferSpectra(MS2ID, idRef$ID_spectra)
        #TODO: put more metrics on the dataframe
        distance <- vapply(seq_along(refSpectra$ptr$id), function(x) {
            a <- .getSpectrum(refSpectra, x)
            row1<- unique(c(a["mass-charge",], massQ))
            row2 <- a["intensity",][match(row1, a["mass-charge",])]
            row1 <- intQ[match(row1, massQ)]
            row1[is.na(row1)] <- 0
            row2[is.na(row2)] <- 0
            lsa::cosine(row1, row2)
        }, FUN.VALUE=3.2)

        hit <- distance >= cosSimThresh
        if(any(hit))
            return(data.frame(REF_idSpectra=refSpectra$ptr$id[hit],
                              distance=distance[hit]))
        else
            return(NA)
    })

    names(distances) <- names(QRY$Spectra)
    #remove query spectra with no hits
    distances <- distances[!is.na(distances)]

    if(length(distances) == 0){
        stop("No query spectra have satisfactory results.")
    }

    #prepare the result
    rslt <- list()
    # crossRef
    rslt$crossRef <-  do.call(rbind, distances)
    rslt$crossRef$QRY_idSpectra <- as.integer(rep(names(distances),
                                       vapply(distances, nrow, FUN.VALUE = 3)))

    SQLwhere <- .appendSQLwhere("ID_spectra",
                                unique(rslt$crossRef$REF_idSpectra), mode="IN")
    crossRef <- .getSQLrecords(MS2ID, "*", "crossRef_SpectrComp",
                                            SQLwhere)
    rslt$crossRef$REF_idCompound <- crossRef$ID_metabolite[
        match(rslt$crossRef$REF_idSpectra, crossRef$ID_spectra)]

    #QRY_metadata
    rslt$QRY_metadata <- QRY$Metadata[match(
        unique(rslt$crossRef$QRY_idSpectra), QRY$Metadata$idSpectra),]
    rslt$QRY_metadata <- dplyr::rename(rslt$QRY_metadata,
                                       QRY_idSpectra=idSpectra)
    #QRY_spectra
    relevantQRYSpectra <- names(QRY$Spectra) %in% unique(rslt$crossRef$QRY_idSpectra)
    rslt$QRY_spectra <- QRY$Spectra
    rslt$QRY_spectra <- rslt$QRY_spectra[relevantQRYSpectra]

    #REF_metaSpectra
    SQLwhere <- .appendSQLwhere("ID_spectra",
                                unique(rslt$crossRef$REF_idSpectra), mode="IN")
    rslt$REF_metaSpectra <- .getSQLrecords(MS2ID, "*", "metaSpectrum",
                                            SQLwhere)
    rslt$REF_metaSpectra <- dplyr::rename(rslt$REF_metaSpectra,
                                           REF_idSpectra=ID_spectra)

    #REF_spectra
    refSpectra <- .bufferSpectra(MS2ID, unique(rslt$crossRef$REF_idSpectra))
    rslt$REF_spectra <- lapply(seq_along(refSpectra$ptr$id), function(x)
        .getSpectrum(refSpectra, x))
    names(rslt$REF_spectra) <- refSpectra$ptr$id

    #REF_metaCompound
    SQLwhere <- .appendSQLwhere("ID_metabolite",
                                unique(rslt$crossRef$REF_idCompound), mode="IN")
    rslt$REF_metaCompound <- .getSQLrecords(MS2ID, "*", "metaCompound",
                                            SQLwhere)
    rslt$REF_metaCompound <- dplyr::rename(rslt$REF_metaCompound,
                                           REF_idCompound=ID_metabolite)

    #obtain possible adducts
    ionizTable <- rslt$QRY_metadata[match(rslt$crossRef$QRY_idSpectra,
                                          rslt$QRY_metadata$QRY_idSpectra),
                                    c("precursorMZ", "polarity")]
    ionizTable$REFMmi <- rslt$REF_metaCompound[match(
        rslt$crossRef$REF_idCompound,
        rslt$REF_metaCompound$REF_idCompound), "REFMmi"]
    rslt$crossRef$propAdduct <- .getAdducts(ionizTable, massError)

    if(cmnNeutralMass){
        cmnNM <- !is.na(rslt$crossRef$propAdduct)
        if(!any(cmnNM)){
            stop("No query spectra have satisfactory results.")
        }
        rslt$crossRef <- rslt$crossRef[cmnNM,]
        rslt$QRY_metadata <- rslt$QRY_metadata[
            rslt$QRY_metadata$QRY_idSpectra %in% rslt$crossRef$QRY_idSpectra,]
        rslt$REF_metaSpectra <- rslt$REF_metaSpectra[
            rslt$REF_metaSpectra$REF_idSpectra %in% rslt$crossRef$REF_idSpectra,]
        rslt$REF_metaCompound <- rslt$REF_metaCompound[
            rslt$REF_metaCompound$REF_idCompound %in%
                rslt$crossRef$REF_idCompound,]
        rslt$QRY_spectra <- rslt$QRY_spectra[
            names(rslt$QRY_spectra) %in% rslt$crossRef$QRY_idSpectra]
        rslt$REF_spectra <- rslt$REF_spectra[
            names(rslt$REF_spectra) %in% rslt$crossRef$REF_idSpectra]
    }

    #obtain number of common masses
    rslt$crossRef$cmnMasses <- vapply(seq_len(nrow(rslt$crossRef)), function(x)
    {
        a <- .getFragments(rslt$QRY_spectra, rslt$crossRef$QRY_idSpectra[x],
                           'mass-charge')
        b <- .getFragments(rslt$REF_spectra, rslt$crossRef$REF_idSpectra[x],
                           'mass-charge')
        sum(a %in% b)
    }, FUN.VALUE = 3)

    runningTime <- round(as.numeric(Sys.time()-annotTime, units = "mins"), 2)
    workVar$annotationDuration <-  paste(runningTime, "min")
    rslt$workVar <- workVar
    return(rslt)
}
#TODO: COMPARE RESULTS WITH PREVOIUS VERSION
#TODO: CHECK MS2ID OBJECT
#TODO: SHINY
#TODO: CONSENSUS SPECTRA
#TODO: INCLUDE MORE METRICS
