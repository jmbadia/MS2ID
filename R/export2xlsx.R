#' @title Export an annotation object to an xlsx file
#'
#' @description \code{export2xlsx} function exports the data contained in an
#' \linkS4class{Annot} object to a xlsx file.
#'
#' @details
#' \code{export2xlsx} function distributes the annotation results among
#' different worksheets of an xlsx file, one per query file, and every row
#' corresponds to an annotation. Annotations of the same query spectrum have
#' identical row colours. Contiguous annotations with equal query precursor mass
#' have similar row colours (e.g. variations of orange); the best distance
#' metric and the most common compound name among those contiguous annotations
#' are bolded in red and black. Moreover:
#' \itemize{\item Column names that contain QRY, REF or CONS refer to query,
#' reference or consensus entities, respectively. e.g. \code{QRYprecursorMz}
#' contains the precursor M/z of the query spectra, \code{REFadduct} the adducts
#' of the reference spectra and \code{QRYrtime_CONS} the retention times of the
#' query spectra that conform every consensus spectrum. \item The column
#' \code{massNum.QRY_CMN_REF} (e.g. 7/4/5) shows the number of fragments of the
#' query and reference spectra, and the number of fragments they have in common.
#' \item \code{REFinchikey} contains the inchikey of the reference compound and
#' a link to a Pubchem's page listing compounds with that inchikey. \item
#' \code{ppmPrecMass} refers to the absolute difference value in ppm between
#' query and reference precursor masses.}
#'
#' @param anRslt \linkS4class{Annot} object with the results to export.
#' @param summarizeHits boolean(1). A TRUE value summarizes the resulting excel,
#'   keeping only, for every query spectrum, the best annotation per compound.
#' @param file char(1) name of the resulting xlsx file
#' @param overwrite boolean(1) If TRUE, overwrite any existing file.
#' @param ... other arguments passed to function
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#' @seealso \linkS4class{Annot} class and \code{\link{annotate}} function.
#' @example man/examples/loadMS2ID.R
#' @example man/examples/selectQuerySpectra.R
#' @examples
#'
#' ## ANNOTATION ---
#' library(MS2ID)
#' MS2IDlib <- MS2ID(MS2IDFolder)
#' annotResult <- annotate(QRYdata = queryFolder, MS2ID = MS2IDlib)
#'
#' ## EXPORT TO XLSX FILE---
#' export2xlsx(anRslt = annotResult, file = "sample",
#'             summarizeHits = FALSE, overwrite=TRUE)
#' @export
export2xlsx <- function(anRslt, file = NULL, summarizeHits = TRUE,
                        overwrite = FALSE, ...){
    argmnts <- c(as.list(environment()), list(...))
    reqClasses <- c(anRslt = "Annot", summarizeHits = "logical",
                    file = "character", overwrite = "logical")
    reqClasses <- reqClasses[names(reqClasses) %in% names(argmnts)]
    .checkTypes(argmnts[match(names(reqClasses), names(argmnts))], reqClasses)
    if(missing(anRslt))
            stop("'anRslt' argument is mandatory")
    dfRslt <- .export2df(anRslt, summarizeHits = summarizeHits, ...)
    #Save xlsx file
    .create_xlsx(data = dfRslt, file = file, overwrite = overwrite, ...)
}

#' @importFrom openxlsx createStyle createWorkbook addWorksheet writeData addStyle writeFormula
.create_xlsx <- function(data = NULL, metric = "cosine", file = NULL, ...){
    argmnts <- c(as.list(environment()), list(...))
    reqClasses <- c(data = "dataframe", metric = "character",
                    file="character")
    reqClasses <- reqClasses[names(reqClasses) %in% names(argmnts)]
    .checkTypes(argmnts[match(names(reqClasses), names(argmnts))], reqClasses)
    if(is.null(data))
        stop("'data' argument is mandatory")
    if(is.null(file))
        stop("'file' argument is mandatory")
    if(length(metric)!=1){
        stop("'metric' must contain ONE string")
    }else if(!(metric %in% c(INCRMETRIC, DECRMETRIC))){
        stop(paste("'metric' must contain ONE the following options:",
                   paste(c(INCRMETRIC, DECRMETRIC), collapse = ", ")))
    }

    #collpse "QRYmassNum","cmnMasses","REFmassNum" in just one column
    massNum_var <- c("QRYmassNum","cmnMasses","REFmassNum")
    data$massNum.QRY_CMN_REF <- vapply(seq_len(nrow(data)), function(nr){
        paste(data[nr, massNum_var], collapse = "/")
    }, FUN.VALUE = "rita")
    data <- dplyr::relocate(data, massNum.QRY_CMN_REF, .before = QRYmassNum)
    data <- dplyr::select(data, !massNum_var)

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
    # m <- which(diff(data$idQRYspect) != 0)

    #Every precursor mass has alternate color
    StyleA1 <- createStyle(border="TopBottom", fgFill = "#ff9a3c")
    StyleA2 <- createStyle(border="TopBottom", fgFill = "#ffc93c")
    StyleB1 <- createStyle(border="TopBottom", fgFill = "#c8f4de")
    StyleB2 <- createStyle(border="TopBottom", fgFill = "#79a8a9")

    lengCol <- seq_len(ncol(data))
    cossimCol <- which(colnames(data) == metric)
    nameCol <- which(colnames(data) == "REFname")
    inchikeyCol <- which(colnames(data) == "REFinchikey")

    ## Create a new workbook
    wb <- createWorkbook(creator = "jmbadia", title="My name here")
    for(idFile in unique(data$QRYdataOrigin)) {
        #select data
        tmp <- data[data$QRYdataOrigin == idFile,]
        #only first 30 characters of the filename
        idFile <- substr(idFile, start = 1, stop = min(30,nchar(idFile)))
        #for every precMass, which has the best cossim
        bst_cossim <- unlist(lapply(unique(tmp$QRYprecursorMz), function(x) {
            maxValue <- max(tmp[tmp$QRYprecursorMz==x, cossimCol])
            which(tmp$QRYprecursorMz==x & tmp[,cossimCol] == maxValue)
        }))
        #for every precMass, which is the most repeated name
        bst_name <- unlist(lapply(unique(tmp$QRYprecursorMz), function(x) {
            mRepN <- names(sort(table(tmp[tmp$QRYprecursorMz == x, "REFname"]),
                              decreasing=TRUE)[1])
            which(tmp$REFname == mRepN & tmp$QRYprecursorMz == x)
        }))

        ## Add a worksheet
        addWorksheet(wb, idFile, gridLines = FALSE)

        ##write data to worksheet 1
        writeData(wb, sheet = idFile, tmp)
        colorA <- unique(tmp$QRYprecursorMz)[c(TRUE, FALSE)]
        colorB <- unique(tmp$QRYprecursorMz)[c(FALSE, TRUE)]

        color1 <- unique(tmp$idQRYspect)[c(TRUE, FALSE)]
        color2 <- unique(tmp$idQRYspect)[c(FALSE, TRUE)]

        clA <- tmp$QRYprecursorMz %in% colorA
        clB <- tmp$QRYprecursorMz %in% colorB
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
        cellInch <- tmp$REFinchikey
        noInch <- is.na(cellInch)
        names(cellInch) <- cellInch
        names(cellInch)[noInch] <- ""
        cellInch[noInch] <- ""
        baselink <- 'https://www.ncbi.nlm.nih.gov/pccompound?term=%22'
        cellInch[!noInch] <- paste0(baselink, cellInch[!noInch],
                                    '%22[InChIKey]')
        class(cellInch) <- "hyperlink"
        openxlsx::writeData(wb, sheet = idFile, x = cellInch, startRow = 2,
                  startCol = inchikeyCol)
    }
    #add xlsx extension if missing
    ex <- strsplit(basename(file), split="\\.")[[1]]
    if(!identical(ex[-1], "xlsx")) file <- paste0(file, ".xlsx")
    openxlsx::saveWorkbook(wb, file = file, ...)
}
