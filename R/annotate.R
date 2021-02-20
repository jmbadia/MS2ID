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
#' @param massErrThresh TODO
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
                     massErrThresh=20, cmnPrecMass=FALSE, cmnNeutralMass=TRUE,
                     cmnPeaks=2, cmnTopPeaks=5, cmnPolarity= TRUE, db="all",
                     nature="all", nsamples, ...){

    #check argument types
    reqClasses <- c(QRYdir="character", MS2ID="MS2ID",
                    noiseThresh="numeric", cmnPeaks="integer",
                    cmnTopPeaks="integer", db="character", nature="character",
                    cmnPolarity= "logical", cmnPrecMass="logical",
                    cmnNeutralMass="logical", cosSimThresh="numeric",
                    massErrThresh="numeric")
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
    if(massErrThresh < 0)
        stop("'massErrThresh' is expected to be a positive number")

    identDate = format(Sys.time(),"%Y%m%d_%H%M")

    message("Loading query spectra ...")
    QRY <- .loadSpectra(dirPath=QRYdir, nsamples=nsamples, ...)
    #QRY <- .loadSpectra(dirPath=QRYdir, nsamples=1000)

    #check queryif argument "cmnPolarity" can be applied
    if(cmnPolarity){
        if(!all(QRY$Metadata$polarity %in% c(0,1)))
            stop("cmnPolarity=TRUE can not be applied because some query spectra
                 have unknown polarity (neither 1 (positive) nor 0 (negative))")
    }
    #remove invalid spectra
    QRY <- .pruneSpectra(QRY)

    # Clean spectra. Remove fragments with intensity < 1% base peak
    # (considering noiseThresh=0.01)
    QRY$Spectra$spectra <- lapply(QRY$Spectra$spectra, function(x) {
        x[,x["intensity",] > noiseThresh * max(x["intensity", ]), drop = F]
    })

    # bin spectra
    #TODO: Make this cte a variable (dependent of the mz index resolution) and
    # check if its consistent with MS2ID bin
    dec2binFrag=2
    message("Binning query spectra ...")
    QRY$Spectra$spectra <- .binSpectra(spectra=QRY$Spectra$spectra,
                                       dec2binFrag)

    # make up a SQL sentence according global restrictions (db, nature,
    SQLwhere <- ""
    if(db!="all")
        SQLwhere <- paste0(SQLwhere, " AND s.REFID_db='", db,"'")
    if(nature!="all")
        SQLwhere <- paste0(SQLwhere, " AND s.REFnature='", nature,"'")

    rslt <- lapply(seq_along(QRY$Spectra$idspctra), function(idQspctr){
        posMetadata <- which(QRY$Metadata$idspctra ==
                                 QRY$Spectra$idspctra[idQspctr])
        mzV <- QRY$Spectra$spectra[[idQspctr]]["mass-charge",]
        idRef <- .getIDref_bymzIndex(mzVector=mzV,  ms2idObj=MS2ID,
                                     cmnPeaks=cmnPeaks,
                                     cmnTopPeaks=cmnTopPeaks)
        idRef <- paste(idRef, collapse = ", ")

        #apply polarity filter
        if(cmnPolarity){
            SQLwhere <- paste0(SQLwhere, " AND s.REFpolarity='",
                               QRY$Metadata$polarity[posMetadata], "'")
        }
        #apply precursor mass filter
        if(cmnPrecMass){
            minPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 - massErrThresh/10^6)
            maxPrecMass <- QRY$Metadata$precursorMZ[posMetadata] *
                (1 + massErrThresh/10^6)
            SQLwhere <- paste0(SQLwhere, " AND s.REFprecursor_mz BETWEEN ",
                               minPrecMass," AND ", maxPrecMass)
        }

        if(!cmnNeutralMass){
            idRef <- DBI::dbGetQuery(MS2ID@dbcon,
                                 paste0("SELECT s.ID_spectra FROM metaSpectrum s
                                       WHERE ID_spectra IN (", idRef,")",
                                       SQLwhere
                                       ))
        }else{
            #TODO: MODIFY ACCORDING cmnNeutralMass requirements
            idRef <- DBI::dbGetQuery(MS2ID@dbcon,
                            paste("SELECT sc.*, s.REFmzRangeEnd, c.REFname
                      FROM crossRef_SpectrComp sc",
                                  "JOIN metaSpectrum s USING(ID_spectra)",
                                  "JOIN metaCompound c USING(ID_metabolite)
                      WHERE s.ID_spectra IN (", idRef,")
                      AND s.REFpolarity=1
                      AND s.REFID_db='Metlin'
                      "))
        }
        #remove Query spectrum with no targeted spectra
       if(nrow(idRef) == 0) return(NULL)

        idRef <- paste(unlist(idRef), collapse = ", ")

        #TODO utilitza function que torni readPos
        spectraPTR <- DBI::dbGetQuery(MS2ID@dbcon,
                               paste("SELECT * FROM spectraPTR
                               WHERE id IN (", idRef,") "))
        #import ALL spectra (on disk -> RAM) to avoid concurrent acces to disk
        readPos <- unlist(lapply(seq_len(nrow(spectraPTR)), function(x)
            seq_len(spectraPTR$numItems[x]) + spectraPTR$startPos[x]))
        refSpectra_thisQ <- MS2ID@spectracon[, readPos, drop=FALSE]

        #recalculate startPos considering refSpectra_thisQ is a subset
        spectraPTR$startPos <- c(0, cumsum(head(spectraPTR$numItems,-1)))

        massm <- QRY$Spectra$spectra[[idQspctr]]["mass-charge",]
        intm <- QRY$Spectra$spectra[[idQspctr]]["intensity",]

        cossim <-vapply(seq_along(spectraPTR$startPos), function(x) {
            a <- refSpectra_thisQ[, seq_along(spectraPTR$numItems[x])+
                                      spectraPTR$startPos[x],
                                    drop=FALSE]
            row1<- unique(c(a["mass-charge",], massm))
            row2 <- a["intensity",][match(row1, a["mass-charge",])]
            row1 <- intm[match(row1, massm)]
            row1[is.na(row1)] <- 0
            row2[is.na(row2)] <- 0
            lsa::cosine(row1, row2)
        }, FUN.VALUE=3.2)

    })
    return(rslt)
}
#TODO: obtain result
#TODO: check all filters
