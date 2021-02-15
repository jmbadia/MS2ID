.getIDref_bymzIndex <- function(mzVector, ms2idObj, cmnPeaks=2, cmnTopPeaks=5){
    mzVector <- paste(mzVector, collapse = ", ")
    #for every mzIndex table, obtain positions of idREF to read
    mzPointers <- lapply(seq_along(MZINDEXPTR), function(i){
        mzIndexPTR <- DBI::dbGetQuery(ms2idObj@dbcon,
                                      paste("SELECT * FROM", MZINDEXPTR[i],
                                            "WHERE id IN (", mzVector,")"))
        unlist(lapply(seq_len(nrow(mzIndexPTR)), function(x)
            seq_len(mzIndexPTR$numItems[x]) + mzIndexPTR$startPos[x]))
    })
    #cmnTopPeaks restriction: Obtain idREF whose first 'cmnTopPeaks'
    # peaks have at list a mzVector
    IDref <- unique(
        ms2idObj@mzIndexcon[unlist(mzPointers[seq_len(cmnTopPeaks)])])
    #cmnPeaks restriction: Keep only idREF with common peaks > cmnPeaks
    IDref <- IDref[vapply(IDref, function(x)
        sum(ms2idObj@mzIndexcon[unlist(mzPointers)] == x) >= cmnPeaks,
        FUN.VALUE=T)]

    return(IDref)
}
