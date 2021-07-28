.mzIdx <- function(DB){
    mzIndex <- list(lstmzIdx = vector("list", 6), spectrIdx = vector("list", 6))
    ## lstmzIdx[[m]] => vector of mz, index to obtain vector of idspectra whose "m" most intense fragment is such mz. m=6 refers to idspectra with such mz but not on its NOT top 5 most intense fragments
    ## lstmzIdx[[2]] == 19.02 19.04 19.10 19.14
    ## spectrIdx[[m]] # list of vectors. Each item is a vector of id spectra whose m most intense fragment is such mz. m=6 refers...
    ## (spectrIdx[[m]][lstmzIdx[[m]] == 19.10] #vector of idspectra whose m most intense fragment is 19.10).
    ## unlist(lapply(seq_len(n), function(x) spectrIdx[[n]][lstmzIdx[[n]] == 19.10]))#vector of idspectra with mz=19.10 on its top n most intense fragments
    # (only) mz of every MS2 spectra sorted by intensity
    spectra_mz <- lapply(DB$fragments$spectra,
                         function(x){
                             xmz <- x["mass-charge",
                                      order(x["intensity",], decreasing = TRUE),
                                      drop = FALSE]
                             return(round(xmz, 2))
                         })
    ## Maxim matrix elemnts affordable
    maxMatrixElemnts <- 3e6

    #loop picking which idspectra has "mz" on its idSub most intense fragment
    for(idSub in seq_along(mzIndex$lstmzIdx)){
        longSpectra <- vapply(spectra_mz, function(x) length(x) >= idSub,
                              FUN.VALUE = TRUE)
        idSpectra_Sub <- DB$fragments$ID_spectra[longSpectra]
        if(idSub != length(mzIndex$lstmzIdx)) { #take most "idSub"nd intense fragment
            spectra_mz_Sub <- lapply(spectra_mz[which(longSpectra)],
                                     function(x) x[idSub])
        } else { #remaining fragments
            spectra_mz_Sub <- lapply(spectra_mz[which(longSpectra)],
                                     function(x) x[-(1:(idSub-1))])
        }
        ## Maxim spectra number feasible per loop
        l1 <- length(unique(unlist(spectra_mz_Sub)))
        l2 <- length(spectra_mz_Sub)
        maxSpectraElemnts <- round(sqrt(maxMatrixElemnts/(l1/l2)))
        spectrIntervals <- c(
            seq(from = 1, to = length(spectra_mz_Sub), by = maxSpectraElemnts),
            1 + length(spectra_mz_Sub))

        ## mz index
        mzIdx <- vector("numeric")
        ## spectra positions (every list corresponds to a mzIdx)
        posSpectrIdx <- vector("list")
        for(iInterval in seq_len(length(spectrIntervals)-1)){
            spIv_1 <- spectrIntervals[iInterval]
            spIv_2 <- spectrIntervals[1 + iInterval]
            mz_Sub2 <- spectra_mz_Sub[spIv_1:(spIv_2 - 1)]
            mzIdx_tmp <- unique(unlist(mz_Sub2)) #unique mz. Points out the row with idspectra
            mm <- matrix(FALSE, nrow = length(mzIdx_tmp),
                         ncol = length(mz_Sub2))
            # positions of mzIdx_tmp corresponding to that mz
            mz_Sub2_pos <- lapply(mz_Sub2, function(x) match(x, mzIdx_tmp))
            for(ia in seq_along(mz_Sub2)){
                mm[mz_Sub2_pos[[ia]], ia ] <- TRUE
            }
            mzIdx <- c(mzIdx, mzIdx_tmp)
            posSpectrIdx <- c(posSpectrIdx,
                              lapply(seq_len(nrow(mm)), function(nr)
                                  spIv_1 - 1 + which(mm[nr,])
                                  )
                              )
        }
        rm(spectra_mz_Sub, mm)
        #summarize duplicated mz results (obtained in different loops)
        mzIdx_unq <- unique(mzIdx)
        posSpectrIdx <- lapply(mzIdx_unq, function(idd)
            unlist(posSpectrIdx[which(mzIdx == idd)])
            ) #15min

        ## sort mzIdx
        posSpectrIdx <- posSpectrIdx[match(sort(mzIdx_unq),mzIdx_unq)]
        mzIndex$lstmzIdx[[idSub]] <- sort(mzIdx_unq)

        ##replace Spectra position with Spectra id
        mzIndex$spectrIdx[[idSub]] <- lapply(posSpectrIdx, function(iPos){
            idSpectra_Sub[iPos]
        }) #15min
        message(paste("Read",idSub,"of",length(mzIndex$lstmzIdx)))
    }
    return(mzIndex)
}
