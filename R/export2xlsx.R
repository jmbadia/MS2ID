#' Export an annotation as xlsx file
#'
#' @param anRslt Annot object with the results to be exporte
#' @param noCmnNeutralMassCos numeric(1) with a cosine threshold to subset Non common neutral mass results.
#' @param file char(1) name (with path) of the exported file
#' @param overwrite boolean(1) If TRUE, overwrite any existing file.
export2xlsx <- function(anRslt, noCmnNeutralMassCos, ...){
    #TODO CHECK ARGUMENTS
    hits <- hits(anRslt)
    if(!missing(noCmnNeutralMassCos)){
        filtVal <- hits$cosine >= cosNOcmnNeutralMass | !is.na(hits$propAdduct)
        hits <- hits[valuableResult,]
    }
    #TODO adapt for others methods than cosine
    decreasMethod <- "cosine" %in% c("topsoe", "fidelity", "squared_chord")
    hits <- hits %>% group_by(idQRYspect, idREFcomp) %>%
        top_n(ifelse(decreasMethod, -1, 1), cosine) %>%
        distinct(idQRYspect, idREFcomp, .keep_all = T) %>%
        arrange(idQRYspect, -cosine)
    REFmetaC <- refCompound(anRslt) %>% rename_with( ~ paste0("REF", .x))
    REFmetaS <- refSpectra(anRslt) %>% Spectra::spectraData() %>%
        as.data.frame() %>% rename_with( ~ paste0("REF", .x))
    REFmetaS$REFmassNum <- vapply(mz(refSpectra(anRslt)), length, FUN.VALUE = 2)
    QRYmetaS <- qrySpectra(anRslt) %>% Spectra::spectraData() %>%
        as.data.frame() %>% rename_with( ~ paste0("QRY", .x))
    QRYmetaS$QRYmassNum <- vapply(mz(qrySpectra(anRslt)), length, FUN.VALUE = 2)
    anRslt <- hits %>%
        merge(REFmetaC, by.x="idREFcomp", by.y="REFid", all.y=F) %>%
        merge(REFmetaS, by.x="idREFspect", by.y="REFid", all.y=F,
              suffixes = c(".comp",".spectra")) %>%
        merge(QRYmetaS, by.x="idQRYspect", by.y="QRYid", all.y=F)

    anRslt <- arrange(anRslt, QRYdataOrigin, QRYrtime, QRYprecursorMZ, -cosine)

    #ORDER & SUBSET columns
    mainVar <- c("QRYprecursorMZ", "QRYrtime", "QRYacquisitionNum", "QRYacquisitionNum_CONS", "REFMmi","propAdduct","REFadduct","REFprecursorMz","cosine","REFname","REFformula","REFinchikey","REFcasNum","QRYmassNum","cmnMasses","REFmassNum","QRYcollisionEnergy","REFCE","QRYpolarity","REFpolarity","QRYprecursorCharge","QRYprecInt","REFnature","REFinstrument","REFinstrumentType", "REFionSource","QRYmsLevel","idQRYspect","idREFspect","idREFcomp","REFID_db.comp","REFID_db.spect","QRYdataOrigin")
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

    #Save xlsx file
    .create_xlsx(data = anRslt, ...)
}

.create_xlsx <- function(data, ...){
    #remove empty columns
    data <- data[ , !vapply(data, function(col){
        all(is.na(col))
    }, FUN.VALUE = T)]

    #concatenate "list" columns
    for(idcol in seq_len(ncol(data))){
        if(is.list(data[,idcol])) {
            data[,idcol] <- vapply(data[,idcol], function(r)
                paste(unlist(r), collapse=", "), FUN.VALUE = "rita")
        }
    }

    ## style for header
    headStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF",
                             halign = "center", fgFill = "#596e79",
                             border="TopBottom", borderColour = "#4F81BD")

    ## style for body
    bodyStyle <- createStyle(border="TopBottom", borderColour = "#4F81BD")
    m <- which(diff(data$idQRYspect)!=0)

    #Every precursor mass has alternate color
    StyleA1 <- createStyle(border="TopBottom", fgFill = "#ff9a3c")
    StyleA2 <- createStyle(border="TopBottom", fgFill = "#ffc93c")
    StyleB1 <- createStyle(border="TopBottom", fgFill = "#c8f4de")
    StyleB2 <- createStyle(border="TopBottom", fgFill = "#79a8a9")

    lengCol <- seq_len(ncol(data))
    cossimCol <- which(colnames(data) %in% c("cosine", "cossim", "distance"))
    nameCol <- which(colnames(data)=="REFname")
    inchikeyCol <- which(colnames(data)=="REFinchikey")

    ## Create a new workbook
    wb <- createWorkbook(creator = "jmbadia", title="My name here")

    noInchikey <- is.na(data$REFinchikey)
    baselink <- 'HYPERLINK("https://www.ncbi.nlm.nih.gov/pccompound?term=%22'
    data$REFinchikey[!noInchikey] <- paste0(baselink,
                                            data$REFinchikey[!noInchikey],
                                            '%22[InChIKey]", "',
                                            data$REFinchikey[!noInchikey],'")')

    for(idFile in unique(data$QRYdataOrigin)) {
        #select data
        tmp <- data[data$QRYdataOrigin==idFile,]

        #only first 30 characters of the filename
        idFile <- substr(idFile, start = 1, stop = min(30,nchar(idFile)))

        #for every precMass, which has the best cossim
        bst_cossim <- unlist(lapply(unique(tmp$QRYprecursorMZ), function(x) {
            maxValue <- max(tmp[tmp$QRYprecursorMZ==x, cossimCol])
            which(tmp$QRYprecursorMZ==x & tmp[,cossimCol] == maxValue)
        }))
        #for every precMass, which is the most repeated name
        bst_name <- unlist(lapply(unique(tmp$QRYprecursorMZ), function(x) {
            mRepN <- names(sort(table(tmp[tmp$QRYprecursorMZ == x, "REFname"]),
                              decreasing=TRUE)[1])
            which(tmp$REFname == mRepN & tmp$QRYprecursorMZ == x)
        }))

        ## Add a worksheet
        addWorksheet(wb, idFile, gridLines = FALSE)

        ##write data to worksheet 1
        writeData(wb, sheet = idFile, tmp)
        #makeHyperlinkString(sheet = idFile, row = 1, col = 1, text = "NULL")
        colorA <- unique(tmp$QRYprecursorMZ)[c(TRUE, FALSE)]
        colorB <- unique(tmp$QRYprecursorMZ)[c(FALSE, TRUE)]

        color1 <- unique(tmp$idQRYspect)[c(TRUE, FALSE)]
        color2 <- unique(tmp$idQRYspect)[c(FALSE, TRUE)]

        clA <- tmp$QRYprecursorMZ %in% colorA
        clB <- tmp$QRYprecursorMZ %in% colorB
        cl1 <- tmp$idQRYspect %in% color1
        cl2 <- tmp$idQRYspect %in% color2

        addStyle(wb, sheet = idFile, headStyle, rows = 1,
                 cols = seq_len(ncol(tmp)), gridExpand = TRUE)
        addStyle(wb, sheet = idFile, StyleA1, rows = 1 + which(clA & cl1),
                 cols = lengCol, gridExpand = TRUE)
        addStyle(wb, sheet = idFile, StyleA2, rows = 1+which(clA & cl2),
                 cols = lengCol, gridExpand = TRUE)
        addStyle(wb, sheet = idFile, StyleB1, rows = 1+which(clB & cl1),
                 cols = lengCol, gridExpand = TRUE)
        addStyle(wb, sheet = idFile, StyleB2, rows = 1+which(clB & cl2),
                 cols = lengCol, gridExpand = TRUE)

        addStyle(wb, sheet = idFile,
                 createStyle(textDecoration = "bold", fontColour = "#d01257"),
                 rows = 1 + bst_cossim, cols = cossimCol, gridExpand = TRUE,
                 stack = TRUE)

        addStyle(wb, sheet = idFile, createStyle(textDecoration = "bold"),
                 rows = 1 + bst_name, cols = nameCol, gridExpand = TRUE,
                 stack = TRUE)

        #add link
        writeFormula(wb, sheet =idFile, startRow = 2, startCol = inchikeyCol,
                     x = tmp$REFinchikey)

    }
    openxlsx::saveWorkbook(wb, ...)
}
