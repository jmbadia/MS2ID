## Utility functions related to adducts analysis

#' Get common numbers considering error
#'
#' @param rawVect numeric(n) ordered and with no NA.
#' @param queryVect numeric(n)
#' @param massError error in ppm to consider
#'
#' @return returns the POSITIONS of rawVect that contain (considering
#' massError) any number of the queryVect.
#' @noRd
.posWhere <- function(rawVect, queryVect, massError=0){
    minVal <- queryVect*(1-massError/10^6)
    maxVal <- queryVect*(1+massError/10^6)
    minPos <- 1+findInterval(minVal, rawVect, left.open=TRUE)
    maxPos <- findInterval(maxVal, rawVect)
    #only intervals where minPos <= maxPos
    contt <- !(minPos > maxPos)
    rslt <- unlist(apply(cbind(minPos[contt],maxPos[contt]), 1,
                         function (x) x[1]:x[2]))
    return(unique(as.vector(rslt)))
}

#' obtain possible adducts
#'
#' @param ionizTable dataframe with Query precursor mass, polarity,
#'  and Reference Mmi columns. Every row describes a ionization
#' @param massError error in ppm to consider
#'
#' @return character(1) with the proposed adducts (collapsed)
#' @noRd
.getAdducts <- function(ionizTable, massError){
    avlblData <- which(!is.na(ionizTable$REFMmi) &
                           !is.na(ionizTable$precursorMZ))
    ionizTable$propAdduct <- NA_character_
    ionizTable$propAdduct[avlblData] <- vapply(avlblData, function(idRow){
        QRYMmi <- .propQMmi(ionizTable$precursorMZ[idRow],
                            ionizTable$polarity[idRow])
        p <- .posWhere(rawVect=sort(QRYMmi),
                       queryVect=ionizTable$REFMmi[idRow],
                       massError=massError)
        if(length(p)!=0) {
            polMatch <- .getPolPos(ionizTable$polarity[idRow])
            adductName <- paste(adducts$Name[polMatch][order(QRYMmi)][p],
                                collapse = " / ")
        }else
            adductName <- NA_character_
        return(adductName)
    }, FUN.VALUE = "julia")
    return(ionizTable$propAdduct)
}

#' Get possible Mmi from a precursor mass
#'
#' @param precMass numeric(1) defining the precursor mass
#' @param polarity integer(1). Optional, negative (0) or positive (1)
#'
#' @return numeric(n) with the proposed Mmi
#' @noRd
.propQMmi <- function(precMass, polarity){
    polMatch <- .getPolPos(polarity)
    QRYMmi <- (precMass-adducts$Mass[polMatch]) *
        abs(adducts$Charge[polMatch])/adducts$Mult[polMatch]
    return(QRYMmi)
}

#return rows of adducts table according polarity
.getPolPos <- function(polarity){
    if(polarity %in% 0:1){
        polMatch <- adducts$Ion_mode == polarity
    }else{
        polMatch <- rep(T, nrow(adducts))
    }
    return(polMatch)
}
