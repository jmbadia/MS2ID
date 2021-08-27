#dec2binFrag = 2 # necessary bin to locate query fragment into the mzIndex
.queryMzIndex <- function(QRYspct, ms2idObj, cmnFrags, dec2binFrag = 2){
    mzVector <- round(
        QRYspct["mass-charge", order(QRYspct['intensity',],
                                     decreasing = TRUE)], dec2binFrag)
    mzVector <- mzVector[seq_len(min(length(mzVector), cmnFrags[2]))]
    #for every mzIndex table, obtain positions of idREF to read
    SQLwhere <- .appendSQLwhere("id", mzVector, mode="IN")
    mzPointe <- lapply(MZINDEXPTR[seq_len(cmnFrags[2])], function(i){
        mzIndexPTR <- .getSQLrecords(ms2idObj, "*", i, SQLwhere)
        unlist(lapply(seq_len(nrow(mzIndexPTR)), function(x)
            seq_len(mzIndexPTR$numItems[x]) + mzIndexPTR$startPos[x]))
    })
    allIDre <- ms2idObj@mzIndexcon[unlist(mzPointe)]
    #cmnFrags[1] restriction: Keep only idREF with common peaks > cmnFrags[1]
    cmnFrags[1] <- min(cmnFrags[1], length(mzVector))
    IDre <- which(tabulate(allIDre) >= cmnFrags[1])
    return(IDre)
}
