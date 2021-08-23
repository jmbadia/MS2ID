#' Cluster spectra to consens
#'
#' Def; spectra group. spectra that has the same main characteristics (mz, polarity, CE, file...).
#' Def; spectra cluster. Spectra (in a group) that is colindant and has a cos>0.8 with the apex spectrum of the cluster
#' Group spectra according its precursor mass, collision energy and polarity.
#' Then subgroup contiguous spectra that have good (cosine) similarity. Each
#' subgroup will form a consensus spectrum.
#'
#' First, all spectra is separated to diffrenet tables with the same Collision
#Energy, polarity, data file and precursor mass. For the last consideration
#first a representative precursor mass is chosen (corresponding to the most
#intense and not already selected precursor); then it is grouped along with all
#its similar precursor mass (considering mass error) For every table we repeat
#the same algorithm. It is selected the most intense precursor 'apex'; its MS2
#spectrum is compared (using cosine similarity method) with the next (forward in
#time) spectrum in time. If they are similar enough it is selected for this
#consensus apex group; if not then the loop stops and the same process is
#repeated but backwards in time. Once done, we apply the loop again to the next
#'apex' found, until no spectrum remains in the table.
#
#' @param s list with spectra and metadata
#' @param consCosThres numeric(1) with  the minimum cosine similarity two contiguous spectra must have to consens a spectrum.
#'
#' @return s with a new column in the metadata dataframe named 'sourceSpect'.
#'   If the row corresponds to: a consensus spectrum, sourceSpect contains a
#'   idSpectra list of the spectra used to consens the spectrum. a primal
#'   spectrum used for a consensus spectrum, sourceSpect is NULL a primal
#'   spectrum non used for a consensus spectrum, sourceSpect is 0
#'
#' @noRd

.cluster2Consens <- function(s, consCosThres = 0.8, massError){
   minScans <- 2 #min number of scans to form a consensus spectra
   # group spectra by mzprecuros
   mdataSplited <- s$Metadata
   mdataSplited$precGroup <- .groupmz(mdataSplited$precursorMZ,
                                      mdataSplited$precursorIntensity,
                                      massError)
   mdataSplited$still <- T #still not avaluated?
   mdataSplited <- group_split(mdataSplited, collisionEnergy, file,
                               precGroup, polarity) #divide in groups
   #cluster spectra
   groupPopul <- vapply(mdataSplited, nrow, FUN.VALUE = 2)
   noGroupSpectra <- vapply(mdataSplited[groupPopul == 1], function(idxG) idxG$idSpectra, FUN.VALUE = 3)
   #for every spectra group with more than 1 spectrum
   spctraClust <- lapply(mdataSplited[groupPopul > 1], function(x){
      genKepp <- vector("list", nrow(x))
      ig <- 1
      while(any(x$still)){
         idAp <- x$idSpectra[x$still][which.max(
            x$precursorIntensity[x$still])[1]]
         apex <- which(x$idSpectra==idAp)#position of the most intense precursor
         keepid <- idAp
         for(i in c(-1, 1)){#backards and then forwards
            cos <- 1
            n <- 0
            while(cos > consCosThres){# check adjacent similarity
               n <- n + 1
               adj <- apex + i *n
               if(adj < 1 | adj > nrow(x)) {
                  break
               }
               if(!x$still[adj]) {
                  break
               }
               idAdj <- x$idSpectra[adj]
               rowdf <- .matchFrag(s$Spectra[[as.character(idAp)]],
                                   s$Spectra[[as.character(idAdj)]],
                                   massError)
               cos <- suppressMessages(
                  philentropy::distance(rowdf, method = "cosine")
                  )
            }
            if(n >= minScans){
               keepid <- c(keepid, x$idSpectra[apex + i * seq_len(n - 1)])
            }
         }
         x$still[x$idSpectra %in% keepid] <- F
         genKepp[[ig]] <- keepid
         ig <- ig+1
      }
      return(genKepp)
   })
   spctraClust <- unlist(spctraClust, recursive = FALSE)
   spctraClust <- spctraClust[!vapply(spctraClust, is.null, FUN.VALUE = T)]
   noClust <- vapply(spctraClust, function(x) length(x) < minScans , FUN.VALUE = T)

   #dataframe with the consensus spectra METADATA
   tmp <- data.frame(idSpectra = seq_len(sum(!noClust)) +
                        max(s$Metadata$idSpectra))
   tmp$sourceSpect <- spctraClust[!noClust]
   refCols <- c("seqNum" , "acquisitionNum" , "spectrumId", "retentionTime")
   refCols <- refCols[refCols %in% names(s$Metadata)]
   for(idRow in seq_len(nrow(tmp))){
      consensdSpectra <- s$Metadata$idSpectra %in% unlist(tmp$sourceSpect[idRow])
      for(idCol in seq_len(ncol(s$Metadata))){
         differentVal <- unique(s$Metadata[consensdSpectra, idCol])
         if(length(differentVal) == 1)
            tmp[idRow, names(s$Metadata)[idCol]] <- differentVal
      }
      tmp$precursorMZ[idRow] <- mean(s$Metadata$precursorMZ[consensdSpectra])
      #collapse reference columns into one
      for(nameCol in refCols){
         newNameCol <- paste0(nameCol,"_CONS")
         tmp[idRow, newNameCol] <- paste(s$Metadata[consensdSpectra, nameCol],
                                      collapse=", ")
      }
   }
   if(nrow(tmp) < 1){
      message("No consensus spectra was formed")
   }else{
      tmp$rol <- 4L
   }
   s$Metadata$rol <- NA_integer_
   s$Metadata$rol[s$Metadata$idSpectra %in% noGroupSpectra] <- 1L
   s$Metadata$rol[s$Metadata$idSpectra %in% unlist(spctraClust[noClust])] <- 2L
   s$Metadata$rol[s$Metadata$idSpectra %in% unlist(tmp$sourceSpect)] <- 3L
   #Add consensus metadata to metadata dataframe
   s$Metadata <- dplyr::bind_rows(tmp, s$Metadata)
   #assign 0 value to spectra used as primal spectra
   s$Metadata$sourceSpect[s$Metadata$idSpectra %in%
                                unlist(spctraClust[noClust])] <- 0
   return(s)
   #rol => 1=spectrum grouped alone, 2=spectra grouped (n>1) but no results in a cluster, 3=spectra grouped(n>1) that results in cluster, 4=cluster spectrum (made of 3 spectra)
}

#' Consens spectra
#'
#' @param s list with spectra and metadata
#' @param massErrorFrag
#' @param minComm numeric(1) minimum ratio of mz presence in order to be present in the final consensus spectrum
#'
#' @return
#' @noRd
.consens <- function(s, massErrorFrag, minComm = 2/3){
   #consensus spectra to be calculated
   toConsens <- s$Metadata$rol == 4L
   spectraRows <- rownames(s$Spectra[[1]])
   consSpectra <- lapply(s$Metadata$sourceSpect[toConsens], function(s2cons){
      #First consens inner mode every spectrum
      spct2cns <- lapply(s2cons, function(x) .getFragments(s$Spectra, x))
      spct2cns <- lapply(spct2cns, function(spct){
         .consensFragm(spct, massErrorFrag = massErrorFrag,
                       mode = "innerSpectra")
      })
      numSpectra <- length(spct2cns)
      #then, consens spectra
      spct2cns <- do.call(cbind, spct2cns)
      spct2cns <- .consensFragm(spct2cns, massErrorFrag = massErrorFrag,
                                minComm = minComm*numSpectra,
                                mode = "interSpectra")
      rownames(spct2cns) <- spectraRows
      mzSort <- sort(spct2cns[1,], index.return = TRUE)$ix
      spct2cns <- spct2cns[, mzSort, drop = FALSE]
      return(spct2cns)
   })
   names(consSpectra) <- s$Metadata$idSpectra[toConsens]
   s$Spectra <- c(consSpectra, s$Spectra)
   return(s)
}

#' Group mz
#'
#' group spectra according its precursor mz +- mass error (from most to less intense)
#'
#' @param mz numeric vector(n) of mz
#' @param int numeric vector(n) with the mz intensities
#' @param massError in ppm
#'
#' @return numeric vector(n) with the group each mz belongs to
#' @noRd
.groupmz <- function(mz, int, massError){
   group <- rep(NA, length(mz)) # group number to where it belongs
   #Find precursorMz representatives
   mErrorMin <- 1 - massError/10^6
   mErrorMax <- 1 + massError/10^6
   n <- 1
   while(any(is.na(group))){
      maxIntPrec <- which.max(int[is.na(group)])[1]
      pM <- mz[is.na(group)][maxIntPrec]
      coPm <- mz < pM*mErrorMax & mz > pM*mErrorMin
      group[coPm] <- n
      n <- n + 1
   }
   return(group)
}

#' match Fragments
#'
#' To compare 2 spectra (i. e. cosine), mz fragments on both spectrum must be considered equal/not equal considering massError
#'
#' @param spectr1 dataframe with mz and intensity as rows, Optionally, if massErrorFrag is missing, spectr1 has also a row with the massError calculated for every mz and a row to alocate spectra 2 intensities
#' @param spectr2 dataframe with mz and intensity as rows
#' @param massErrorFrag in ppm.
#'
#' @return dataframe with 2 rows, spectr1 & spectr2 intensities. Both in the same column when the mz has been considered equal
#' @noRd
.matchFrag <- function(spectr1, spectr2, massErrorFrag){
   if(!missing(massErrorFrag)){
      spectr1 <- rbind(spectr1,
                       error = spectr1["mass-charge",] * massErrorFrag * 1e-6,
                       intSpectr2 = 0)
   }else{
      spectr1["intSpectr2", ] <- 0
   }
   for(idxrmz in seq_len(ncol(spectr2))){
      mzfragm <- spectr2["mass-charge", idxrmz]
      qnear <- which.min(abs(spectr1["mass-charge", ] - mzfragm))
      near <- spectr1["error", qnear] >= abs(spectr1["mass-charge",
                                                     qnear] - mzfragm)
      if(near){
         spectr1["intSpectr2", qnear] <- spectr1["intSpectr2", qnear] +
            spectr2["intensity", idxrmz]
      }else{
         spectr2["mass-charge", idxrmz] <- NA
      }
   }
   isnaSpect2 <- is.na(spectr2["mass-charge",])
   struct <- cbind(spectr1[c("intensity","intSpectr2"),],
                   rbind(rep(0, sum(isnaSpect2)),
                         spectr2["intensity", isnaSpect2]))
   return(struct)
}

#' match mz
#'
#' Considering two mz vectors (mz1 and mz2) and a mass error, obtain which mz1 equals every mz2 value.
#'
#' @param mz1 numeric(n) vector
#' @param mz1 numeric(n) vector
#' @param massError in ppm.
#'
#' @return numeric(n) vector with the position of mz1 that equals every mzm value (similar to the match() return). E.g. .matchMz(c(1, 7), c(7, 9, 11))=c(2, NA, NA)
#' @noRd
.matchMz <- function(mz1, mz2, massError){
   mz1_err <- mz1 * massError * 1e-6
   for(idxrmz in seq_along(mz2)){
      mzfragm <- mz2[idxrmz]
      qnear <- which.min(abs(mz1 - mzfragm))
      near <- mz1_err[qnear] >= abs(mz1[qnear] - mzfragm)
      mz2[idxrmz] <- ifelse(near, qnear, NA)
   }
   return(mz2)
}

#' consens a spectrum or spectra
#'
#' @param spct spectrum or spectra (merged) to be consensued
#' @param massErrorFrag
#' @param minComm integer(1) rfering to minimum iterations of a mz in primal spectra in order to be present in the cosensus spectrum. Only considered in interSpectra mode
#' @param mode 'innerSpectra' or 'interSpectra'. The former is to bin (sum up) close mz in a spectrum. The last finds common mz inter spectra.
#'
#' @return
#' @noRd
.consensFragm <- function(spct, massErrorFrag, minComm, mode){
   groups <- .groupmz(spct[1,], spct[2,], massErrorFrag)
   ngroups <- vapply(groups, function(x) sum(groups==x), FUN.VALUE = 3)
   if(mode == "innerSpectra"){
      grouped <- ngroups > 1
      if(any(grouped)){
         a <- vapply(unique(groups[grouped]), function(gr){
            thisgroup <- groups == gr
            mz_v <- spct[1 , thisgroup]
            i_v <- spct[2 , thisgroup]
            mz_gr <- sum(i_v * mz_v)/sum(i_v)#weighted mz mean
            i_gr <- sum(i_v)#sum int
            return(c(mz_gr, i_gr))
         }, FUN.VALUE = c(1.2, 1.2))
         spct <- cbind(a, spct[, !grouped, drop=F])
      }
   }else if(mode == "interSpectra"){
      grouped <- ngroups >= minComm
      if(any(grouped)){
         spct <- vapply(unique(groups[grouped]), function(gr){
            thisgroup <- groups == gr
            mz_v <- spct[1 , thisgroup]
            i_v <- spct[2 , thisgroup]
            mz_gr <- sum(i_v * mz_v)/sum(i_v)#weighted mz mean
            i_gr <- sum(i_v)/sum(thisgroup)#mean
            return(c(mz_gr, i_gr))
         }, FUN.VALUE = c(1.2, 1.2))
      }else{
         spct <- NULL
      }
   }
   return(spct)
}

# Considering a dataframe with columns named x and x'_CONS', copies the non NA rows of x'_CONS' into x and removes the x_'CONS' columns. OPtionally resizes  the size of the x'_CONS' content
.mergeCONS <- function(dfCONS, fontsize){
   if(nrow(dfCONS) > 0 & any(grepl("_CONS$", names(dfCONS)))){
      am <- dfCONS %>% select(ends_with("_CONS")) %>%
         rename_with(~ gsub("_CONS", "", .x, fixed = TRUE))
      nuclearname <- names(am)
      CONSna <- is.na(am[, 1])
      if(!missing(fontsize)){
         for(idxcol in seq_len(ncol(am))){
            am[!CONSna, idxcol] <- paste0(
               "<font size=", fontsize, ">", am[!CONSna, idxcol],"</font>")
         }
      }
      dfCONS <- mutate(dfCONS, across(nuclearname, as.character))
      dfCONS[!CONSna, match(nuclearname, names(dfCONS))] <- am[!CONSna,]
      dfCONS <- select(dfCONS, !ends_with("_CONS"))
   }
   return(dfCONS)
}

