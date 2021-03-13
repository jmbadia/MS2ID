#' Structure the annotation results
#'
#' On the basis of the annotation results, the '.processAnnotation()' function
#'  adds metadata, calculates proposed adducts and common masses for every
#'  hit (i.e. every pair query - reference spectrum with satisfactory
#'  distance metrics) and packs altogether into an AnnotdSpectra object.
#'
#' @param dist annotation results
#' @param ms2id ms2id object with the reference library used
#' @param qry query spectra
#' @param mError Mass error specified by the user
#' @param cmnNtMass boolean indicating if hits must be subsetted to those with
#'  proposed adducts explaining the hit
#' @param workVar arguments specified by the user in the annotate function.
#'
#' @return an AnnotdSpectra object with the results of the annotation
#' @noRd
.processAnnotation<- function(dist, ms2id, qry, mError, cmnNtMass, workVar) {
    hits <- do.call(rbind, dist)
    hits$idQRYspect <- as.numeric(rep(names(dist),
                                      vapply(dist, nrow, FUN.VALUE = 3)))

    #crossRef (REF spectra vs REF compound)
    SQLwhere <- .appendSQLwhere("ID_spectra", unique(hits$idREFspect),
                                mode="IN")
    crossRef <- .getSQLrecords(ms2id, "*", "crossRef_SpectrComp", SQLwhere)
    hits$idREFcomp <- crossRef$ID_metabolite[match(hits$idREFspect,
                                                   crossRef$ID_spectra)]
    #compound Metadata
    SQLwhere <- .appendSQLwhere("ID_metabolite", unique(crossRef$ID_metabolite),
                                mode="IN")
    REFcomp <- .getSQLrecords(ms2id, "*", "metaCompound", SQLwhere)
    REFcomp <- dplyr::rename_with(REFcomp, ~ gsub("REF", "", .x, fixed = TRUE))
    REFcomp <- dplyr::rename(REFcomp, id = 'ID_metabolite')

    #obtain possible adducts
    ionizTable <- qry$Metadata[match(hits$idQRYspect, qry$Metadata$idSpectra),
                               c("precursorMZ", "polarity")]
    ionizTable$Mmi <- REFcomp[match(hits$idREFcomp, REFcomp$id), "Mmi"]
    hits$propAdduct <- .getAdducts(ionizTable, mError)

    if(cmnNtMass){
        cmnNM <- !is.na(hits$propAdduct)
        if(!any(cmnNM))  stop("No query spectra have satisfactory results.")
        hits <- hits[cmnNM,]
    }

    #QRYspect
    QRYspect <- qry$Metadata[qry$Metadata$idSpectra %in% hits$idQRYspect,]
    QRYspect <- dplyr::rename(QRYspect, id = 'idSpectra', dataOrigin = 'file',
                              rtime = 'retentionTime')
    orderSpectra <- match(QRYspect$id, names(qry$Spectra))

    QRYspect$mz <- lapply(orderSpectra, function(x)
        qry$Spectra[[x]]['mass-charge',])
    QRYspect$intensity <- lapply(orderSpectra, function(x)
        qry$Spectra[[x]]['intensity',])

    #REFspect
    SQLwhere <- .appendSQLwhere("ID_spectra", unique(hits$idREFspect),
                                mode="IN")
    REFspect <- .getSQLrecords(ms2id, "*", "metaSpectrum", SQLwhere)
    REFspect <- dplyr::rename_with(REFspect,
                                   ~ gsub("REF", "", .x, fixed = TRUE))
    REFspect$msLevel <- as.integer(gsub("MS", "", REFspect$msLevel,
                                        fixed = TRUE))
    REFspect <- dplyr::rename(REFspect, id = 'ID_spectra',
                              precursorMz = 'precursor_mz')
    refSpectra <- .bufferSpectra(ms2id, unique(hits$idREFspect))
    REFspect$mz <- lapply(REFspect$id, function(x){
        .getSpectrum(refSpectra, which(refSpectra$ptr$id==x))['mass-charge',]
    })
    REFspect$intensity <- lapply(REFspect$id, function(x){
        .getSpectrum(refSpectra, which(refSpectra$ptr$id==x))['intensity',]
    })

    #REFcomp
    REFcomp <- REFcomp[REFcomp$id %in% hits$idREFcomp, ]

    #obtain number of common masses
    hits$cmnMasses <- vapply(seq_len(nrow(hits)), function(x) {
        a <- QRYspect$mz[match(hits$idQRYspect[x], QRYspect$id)]
        b <- REFspect$mz[match(hits$idREFspect[x], REFspect$id)]
        sum(unlist(a) %in% unlist(b))
    }, FUN.VALUE = 3)

    REFspect <- Spectra::Spectra(REFspect)
    QRYspect <- Spectra::Spectra(QRYspect)

    runningTime <- round(as.numeric(Sys.time() - workVar$annotationTime,
                                    units = "mins"), 2)
    workVar$annotationDuration <-  paste(runningTime, "min")
    workVar$annotationTime <- as.character(workVar$annotationTime)

    rslt <- AnnotdSpectra(qrySpectra = QRYspect, refSpectra = REFspect,
                          refCompound = REFcomp, hits = hits,
                          infoAnnotation = workVar)
    return(rslt)
}
