#' Extract metadata from an annotation object in dataframe format
#'
#' @param anRslt Annot object with the results to be exported
#' @param metric,metricThresh char(1) and numeric(1) with the metric (and its
#'   threshold respectively) that subsets annotations on the final excel file.
#'   \code{metric} parameter is also used to sort the annottations.
#' @param summarizeHits boolean(1) if TRUE, resulting dataframe keeps only the best hits of every query spectra - Compound pair.
#' @param file char(1) name (with path) of the exported file
#' @param overwrite boolean(1) If TRUE, overwrite any existing file.
#' @noRd
#'
#' @import dplyr
.export2df <- function(anRslt, metric = "cosine", metricThresh,
                       summarizeHits = TRUE){
    #CHECK ARGUMENTS
    argmnts <- c(as.list(environment()))
    reqClasses <- c(anRslt = "Annot", summarizeHits = "logical",
                    metric = "character", metricThresh = "numeric")
    reqClasses <- reqClasses[names(reqClasses) %in% names(argmnts)]
    .checkTypes(argmnts[match(names(reqClasses), names(argmnts))], reqClasses)
    if(missing(anRslt))
        stop("'anRslt' argument is mandatory")
    if(length(metric)!=1){
        stop("'metric' must contain ONE string")
    }else if(!(metric %in% c(INCRMETRIC, DECRMETRIC))){
        stop(paste("'metric' must contain ONE the following options:",
                   paste(c(INCRMETRIC, DECRMETRIC), collapse = ", ")))
    }
    hits <- hits(anRslt)
    if(!(metric %in% names(hits)))
       stop(paste0("Metric '", metric, "' not found in the annotation object"))

    if(!missing(metricThresh)){
        filtVal <- hits[, metric] >= metricThresh
        hits <- hits[filtVal, ]
    }

    #summarize if asked for
    sgnMetric <- ifelse(metric %in% DECRMETRIC, 1, -1)
    if(summarizeHits){
        hits <- hits %>%
            group_by(idQRYspect, idREFcomp) %>%
            top_n(-sgnMetric, !!as.name(metric)) %>%
            distinct(idQRYspect, idREFcomp, .keep_all = T) %>%
            ungroup()
    }

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
        merge(QRYmetaS, by.x="idQRYspect", by.y="QRYid", all.y=F) %>%
        mutate(REFpredicted = replace(REFpredicted, REFpredicted == 0,
                                      "experimental"),
               REFpredicted = replace(REFpredicted, REFpredicted == 1,
                                      "insilico")
               )
    ppmPrecMass <- (anRslt$QRYprecursorMz - anRslt$REFprecursorMz) /
        anRslt$QRYprecursorMz * 1e6
    if(!all(is.na(ppmPrecMass)))
        anRslt$ppmPrecMass <- round(ppmPrecMass, 1)
    #sort columns
    mainVar <- c("QRYprecursorMz","idQRYspect", "idREFspect", "idREFcomp",
                 INCRMETRIC, DECRMETRIC, "QRYmassNum", "cmnMasses",
                 "REFmassNum", "propAdduct", "REFname","REFformula",
                 "REFinchikey", "ppmPrecMass", "QRYcollisionEnergy_txt",
                 "REFcollisionEnergy", "REFexactmass","REFprecursorMz",
                 "REFadduct","REFpredicted", "REFinstrument",
                 "REFinstrumentType","REFID_db.comp", "REFID_db.spectra",
                 "QRYrtime", "QRYacquisitionNum", "QRYacquisitionNum_CONS",
                 "REFCASRN","QRYpolarity", "REFpolarity",
                 "QRYprecursorCharge","QRYprecursorIntensity",
                 "REFionSource", "QRYmsLevel",  "QRYdataOrigin")
    anRslt <- select(anRslt, c(mainVar[mainVar %in% names(anRslt)],
                               names(anRslt)[!names(anRslt) %in% mainVar]))
    #sort rows
    if(summarizeHits){
        anRslt <- anRslt %>%
            arrange(QRYdataOrigin, QRYprecursorMz, idQRYspect,
                    sgnMetric*!!as.name(metric))
    }else{
        anRslt <- anRslt %>%
            group_by(idQRYspect, idREFcomp) %>%
            mutate(top = ifelse(metric %in% DECRMETRIC,
                                min(!!as.name(metric)),
                                max(!!as.name(metric)))) %>%
            arrange(QRYdataOrigin, QRYprecursorMz, idQRYspect,
                    sgnMetric*top, idREFcomp,
                    sgnMetric*!!as.name(metric)) %>%
            ungroup() %>%
            select(-top) %>%
            as.data.frame()
    }
    return(anRslt)
}
