#' Extract metadata from an annotation object in dataframe format
#'
#' @param anRslt Annot object with the results to be exporte
#' @param noCmnNeutralMassDist numeric(1) with a distance threshold to subset Non common neutral mass results.
#' @param metric char(1) with the name of the metric that must subset (according to the distance threshold) and order the identifications.
#' @param file char(1) name (with path) of the exported file
#' @param overwrite boolean(1) If TRUE, overwrite any existing file.
#' @noRd
#'
.export2df <- function(anRslt, noCmnNeutralMassDist, metric = "cosine", ...){
    #TODO CHECK ARGUMENTS

    if(length(metric)!=1){
        stop("'metric' must contain ONE string")
    }else if(!(metric %in% c(INCRMETRIC, DECRMETRIC))){
        stop(paste("'metric' must contain ONE the following options:",
                   paste(c(INCRMETRIC, DECRMETRIC), collapse = ", ")))
    }

    hits <- hits(anRslt)
    if(!(metric %in% names(hits)))
       stop(paste0("Metric '", metric, "' not found in the annotation object"))

    if(!missing(noCmnNeutralMassDist)){
        filtVal <- hits[, metric] >= noCmnNeutralMassDist |
            !is.na(hits$propAdduct)
        hits <- hits[valuableResult,]
    }

    sgnMetric <- ifelse(metric %in% DECRMETRIC, -1, 1)
    hits <- hits %>% dplyr::group_by(idQRYspect, idREFcomp) %>%
        dplyr::top_n(sgnMetric, !!as.name(metric)) %>%
        dplyr::distinct(idQRYspect, idREFcomp, .keep_all = T) %>%
        dplyr::arrange(idQRYspect, sgnMetric*!!as.name(metric))
    REFmetaC <- refCompound(anRslt) %>% rename_with( ~ paste0("REF", .x))
    REFmetaS <- refSpectra(anRslt) %>% Spectra::spectraData() %>%
        as.data.frame() %>% rename_with( ~ paste0("REF", .x))
    REFmetaS$REFmassNum <- vapply(Spectra::mz(refSpectra(anRslt)),
                                  length, FUN.VALUE = 2)
    QRYmetaS <- qrySpectra(anRslt) %>% Spectra::spectraData() %>%
        as.data.frame() %>% rename_with( ~ paste0("QRY", .x))
    QRYmetaS$QRYmassNum <- vapply(Spectra::mz(qrySpectra(anRslt)),
                                  length, FUN.VALUE = 2)
    anRslt <- hits %>%
        merge(REFmetaC, by.x="idREFcomp", by.y="REFid", all.y=F) %>%
        merge(REFmetaS, by.x="idREFspect", by.y="REFid", all.y=F,
              suffixes = c(".comp",".spectra")) %>%
        merge(QRYmetaS, by.x="idQRYspect", by.y="QRYid", all.y=F)

    anRslt <- arrange(anRslt, QRYdataOrigin, QRYrtime, QRYprecursorMZ,
                      sgnMetric*!!as.name(metric))

    #ORDER & SUBSET columns
    mainVar <- c("QRYprecursorMZ", "QRYrtime", "QRYacquisitionNum",
        "QRYacquisitionNum_CONS", "REFMmi","propAdduct","REFadduct",
        "REFprecursorMz", INCRMETRIC, DECRMETRIC, "REFname","REFformula",
        "REFinchikey", "REFcasNum","QRYmassNum","cmnMasses","REFmassNum",
        "QRYcollisionEnergy", "REFCE","QRYpolarity","REFpolarity",
        "QRYprecursorCharge","QRYprecInt", "REFnature","REFinstrument",
        "REFinstrumentType", "REFionSource", "QRYmsLevel","idQRYspect",
        "idREFspect","idREFcomp","REFID_db.comp", "REFID_db.spect",
        "QRYdataOrigin")
    anRslt <- dplyr::select(anRslt,
                            c(mainVar[mainVar %in% names(anRslt)],
                              names(anRslt)[!names(anRslt) %in% mainVar]))

    #collpse "QRYmassNum","cmnMasses","REFmassNum" in just one column
    massNum_var <- c("QRYmassNum","cmnMasses","REFmassNum")
    anRslt$QRY_CMN_REF_massNum <- vapply(seq_len(nrow(anRslt)), function(nr){
        paste(anRslt[nr, massNum_var], collapse = "/")
    }, FUN.VALUE = "rita")
    anRslt <- dplyr::relocate(anRslt, QRY_CMN_REF_massNum, .before=QRYmassNum)
    anRslt <- dplyr::select(anRslt, !massNum_var)
    return(anRslt)
}
