#' Structure the annotation results
#'
#' On the basis of the annotation results, the '.processAnnotation()' function
#'  adds metadata, calculates proposed adducts and common masses for every
#'  hit (i.e. every pair query - reference spectrum with satisfactory
#'  distance metrics) and packs altogether into an Annot object.
#'
#' @param dist annotation results
#' @param ms2id ms2id object with the reference library used
#' @param qry query spectra
#' @param massErrMs1 numeric(1). Mass error to consider in operations with first spectrometer measures, e.g. infere adducts from precursor and natural mass
#' @param massErrMs1 numeric(1). Mass error to consider in operations with non first spectrometer measures (typycally MS2), e.g. calculate common fragments in query - reference spectra pairs
#' @param cmnNtMass boolean indicating if hits must be subsetted to those with
#'  proposed adducts explaining the hit
#' @param workVar arguments specified by the user in the annotate function.
#'
#' @return an Annot object with the results of the annotation
#' @noRd
.processAnnotation<- function(dist, ms2id, qry, lft, massErrMs1, massErrMsn,
                              cmnNtMass, workVar) {
    hits <- do.call(rbind, dist)
    hits$idQRYspect <- as.numeric(rep(names(dist),
                                      vapply(dist, nrow, FUN.VALUE = 3)))
    #crossRef (REF spectra vs REF compound)
    SQLwhere <- .appendSQLwhere("ID_spectra", unique(hits$idREFspect),
                                mode="IN")
    crossRef <- .getSQLrecords(ms2id, "*", "crossRef_SpectrComp", SQLwhere)

    # Check if idREFspect is related to more than one metabolite. In that
    # case replicate n times the row with that idREFspect in order to assign
    #  a idREFcomp to each row
    cmpndSpctrum <- lapply(hits$idREFspect, function(x)
        crossRef[crossRef$ID_spectra == x, "ID_compound"])
    #copy n times the row where n is the number of metabolites per REF spectrum
    hits <- hits[rep(seq_along(cmpndSpctrum),
                      vapply(cmpndSpctrum, length, FUN.VALUE = 1)),]
    hits$idREFcomp <- unlist(cmpndSpctrum)

    #compound Metadata
    SQLwhere <- .appendSQLwhere("ID_compound", unique(crossRef$ID_compound),
                                mode="IN")
    REFcomp <- .getSQLrecords(ms2id, "*", "metaCompound", SQLwhere)
    REFcomp <- dplyr::rename_with(REFcomp, ~ gsub("REF", "", .x, fixed = TRUE))
    REFcomp <- dplyr::rename(REFcomp, id = 'ID_compound')

    #obtain possible adducts
    ionizTable <- qry$Metadata[match(hits$idQRYspect, qry$Metadata$idSpectra),
                               c("precursorMZ", "polarity")]
    ionizTable$Mmi <- REFcomp[match(hits$idREFcomp, REFcomp$id), "exactmass"]
    hits$propAdduct <- .getAdducts(ionizTable, massErrMs1)

    if(cmnNtMass){
        cmnNM <- !is.na(hits$propAdduct)
        if(!any(cmnNM))
            stop("No query spectrum meets the cmnNtMass condition")
        hits <- hits[cmnNM,]
    }
    #QRYspect
    idCONShits <- unique(hits$idQRYspect)
    idSRCShits <- qry$Metadata[qry$Metadata$idSpectra %in% idCONShits,
                               "sourceSpect"]
    idSRCShits <- unique(unlist(idSRCShits))
    QRYspect <- qry$Metadata[qry$Metadata$idSpectra %in% idCONShits, ]
    orderSpectra <- match(as.character(QRYspect$id), names(qry$Spectra))
    QRYspect$mz <- lapply(orderSpectra, function(x)
        qry$Spectra[[x]]['mass-charge',])
    QRYspect$intensity <- lapply(orderSpectra, function(x)
        qry$Spectra[[x]]['intensity',])
    if(!anyNA(lft, recursive=FALSE)){
        LFTspect <- lft$Metadata[lft$Metadata$idSpectra %in% idSRCShits, ]
        orderSpectra <- match(LFTspect$id, names(lft$Spectra))

        LFTspect$mz <- lapply(orderSpectra, function(x)
            lft$Spectra[[x]]['mass-charge',])
        LFTspect$intensity <- lapply(orderSpectra, function(x)
            lft$Spectra[[x]]['intensity',])

        QRYspect <- rbind(QRYspect, LFTspect)
    }
    # old_id = 'id',
    lookup <- c(idOld = 'id', id = 'idSpectra', dataOrigin = 'file',
                precScanNum = "precursorScanNum", precursorMz ="precursorMZ",
                isolationWindowTargetMz = "isolationWindowTargetMZ")
    QRYspect <- QRYspect %>%
        dplyr::rename_with(~ gsub("retentionTime", "rtime", .x)) %>%
        dplyr::rename(any_of(lookup))
    #REFspect
    SQLwhere <- .appendSQLwhere("ID_spectra", unique(hits$idREFspect),
                                mode="IN")
    REFspect <- .getSQLrecords(ms2id, "*", "metaSpectrum", SQLwhere)
    REFspect <- dplyr::rename_with(REFspect,
                                   ~ gsub("REF", "", .x, fixed = TRUE))
    REFspect <- dplyr::rename(REFspect, id = 'ID_spectra')
    refSpectra <- .bufferSpectra(ms2id, unique(hits$idREFspect))
    REFspect$mz <- lapply(REFspect$id, function(x){
        .getSpectrum(refSpectra, which(refSpectra$ptr$id==x))['mass-charge',]
    })
    REFspect$intensity <- lapply(REFspect$id, function(x){
        .getSpectrum(refSpectra, which(refSpectra$ptr$id==x))['intensity',]
    })

    #REFcomp
    REFcomp <- REFcomp[REFcomp$id %in% hits$idREFcomp, ]
    #obtain which query fragments corresponds to each reference fragment
    hits$cmnMasses <- vapply(seq_len(nrow(hits)), function(x) {
        tmpqry <- QRYspect$mz[match(hits$idQRYspect[x], QRYspect$id)] %>%
            unlist()
        tmpref <- REFspect$mz[match(hits$idREFspect[x], REFspect$id)] %>%
            unlist()
        cmnM <- sum(!is.na(.matchMz(tmpqry, tmpref, massErrMsn)))
        return(cmnM)
    }, FUN.VALUE = 3)

    #remove NA columns
    QRYspect <- QRYspect[, !vapply(QRYspect, function(col) all(is.na(col)),
                              FUN.VALUE = T)]
    REFspect <- REFspect[, !vapply(REFspect, function(col) all(is.na(col)),
                                   FUN.VALUE = T)]

    #Spectra package only admits numeric CollisionEnergy values
    if("collisionEnergy" %in% names(QRYspect))
        QRYspect$collisionEnergy <- as.numeric(QRYspect$collisionEnergy)
    else
        QRYspect$collisionEnergy <- NA_real_

    if("collisionEnergy" %in% names(REFspect)){
        REFspect$collisionEnergy_txt <- REFspect$collisionEnergy
        REFspect$collisionEnergy <- suppressWarnings(
            as.numeric(REFspect$collisionEnergy))
    }else
        REFspect$collisionEnergy_txt <- REFspect$collisionEnergy <- NA_real_

    #Convert to spectra object
    QRYspect <- Spectra::Spectra(QRYspect)
    REFspect <- Spectra::Spectra(REFspect)

    runningTime <- round(as.numeric(Sys.time() - workVar$annotationTime,
                                    units = "mins"), 2)
    workVar$annotationDuration <-  paste(runningTime, "min")
    workVar$annotationTime <- as.character(workVar$annotationTime)

    rslt <- Annot(qrySpectra = QRYspect, refSpectra = REFspect,
                  refCompound = REFcomp, hits = hits, infoAnnotation = workVar)
    return(rslt)
}
