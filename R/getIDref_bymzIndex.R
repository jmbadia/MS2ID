.queryMzIndex <- function(mzVector, ms2idObj, cmnPeaks=2, cmnTopPeaks=5){
    #for every mzIndex table, obtain positions of idREF to read
    SQLwhere <- .appendSQLwhere("id", mzVector, mode="IN")
    mzPointers <- lapply(seq_along(MZINDEXPTR), function(i){
        mzIndexPTR <- .getSQLrecords(ms2idObj, "*", MZINDEXPTR[i], SQLwhere)
        unlist(lapply(seq_len(nrow(mzIndexPTR)), function(x)
            seq_len(mzIndexPTR$numItems[x]) + mzIndexPTR$startPos[x]))
    })
    #cmnTopPeaks restriction: Obtain idREF whose first 'cmnTopPeaks'
    # peaks have at list a mzVector
    IDref <- unique(
        ms2idObj@mzIndexcon[unlist(mzPointers[seq_len(cmnTopPeaks)])])
    allIDref <- ms2idObj@mzIndexcon[unlist(mzPointers)]
    #cmnPeaks restriction: Keep only idREF with common peaks > cmnPeaks
    IDref <- which(tabulate(allIDref[allIDref %in% IDref]) >= cmnPeaks)
    return(IDref)
}
