#' Adapt the metadata to the MS2ID structure
#'
#' @param metadata dataframe containing spectra & metabolite metadata (variable
#'   = columns). No id variables must be included (ID_spectra, ID_compound,
#'   ID_db)
#' @param fragm list where each item is a MS2 spectra. Row order must match with
#'   metadata row order
#' @param varsToParse dataframe with the MS2ID variable names as rownames, the
#'   BD original names and the type. Must match with the columns of the metadata
#' @param nameDB char(1) with the name we want for the DB we are parsing
#' @param numCmp int(1) number of compounds in the original DB
#' @param numMs int(1) number of spectra in the original DB
#' @param removeRedCompounds boolean(1) eliminate redundant compounds based
#' on the inchikeys
#' @noRd
.basicMS2IDstrct <- function(metadata, fragm, varsToParse,
                             nameDB="unknown", numCmp, numMs,
                             removeRedCompounds = FALSE)
  {
  #rename variables with MS2ID names
  MS2IDname <- match(colnames(metadata), varsToParse$originalNames)
  colnames(metadata)[!is.na(MS2IDname)] <-
    varsToParse$MS2IDname[MS2IDname[!is.na(MS2IDname)]]
  included <- colnames(metadata)[!is.na(MS2IDname)]
  message(glue::glue("
    \nInput metadata considered:\n \\
    {glue::glue_collapse(included, ', ', last = ' and ')}"))
  if(any(is.na(MS2IDname)))
  {
    notIncluded <- colnames(metadata)[is.na(MS2IDname)]
    message(glue::glue("
    \nInput metadata discarded:\n \\
    {glue::glue_collapse(notIncluded, ', ', last = ' and ')}\\
    (variables not covered by the MS2ID class)"))
  }
  if(!all(varsToParse$MS2IDname %in% varsToParse$originalNames)){
    notIncluded <- varsToParse$MS2IDname[
      !varsToParse$MS2IDname %in% colnames(metadata)]
    message(glue::glue("
     \nInput metadata not present:\n \\
     {glue::glue_collapse(notIncluded, ', ', last = ' and ')}\\
     (variables covered by the MS2ID class but not present in the\\
     input metadata)"))
  }

  #remove not considered variables (even spectrum_id, which we do not need anymore)
  metadata <- metadata[,!is.na(MS2IDname)]
  #remove info from variables we will not use
  varsToParse <- varsToParse[MS2IDname[!is.na(MS2IDname)],]

  #Check data types--------
  classOriginal <- sapply(seq_len(ncol(metadata)), function(idCol)
    class(metadata[,idCol]))
  classMS2ID <- varsToParse$MS2IDdataType[match(colnames(metadata),
                                          varsToParse$MS2IDname)]
  if(!identical(classOriginal, classMS2ID)){
    stop("The following variables do not have the proper data type and must be converted: ")
    notEqualClass <- which(classOriginal!=classMS2ID)
    for(idC in notEqualClass){
      print(paste(colnames(metadata)[idC], ":",
                  classOriginal[idC], "=>", classMS2ID[idC]))
    }
  }

  #0. Check & repair----------------
  message("\nChecking & repairing data------")
#remove entries with NA spectra
NASpectra <- is.na(fragm)
metadata <- metadata[!NASpectra,]
fragm <- fragm[!NASpectra]
#assign NA to empty entries
metadata[metadata == ""] <- NA_character_
# Count NA / variable
na_count <- sapply(metadata, function(y) sum(length(which(is.na(y)))))
showNA <- paste(round(100*na_count/nrow(metadata)),"%")# % each variable's NA)
names(showNA) <- names(na_count)
message(paste("NA presence on every variable:\n", .print_and_capture(showNA)))

message("\n-------- 2/n. REDUNDANT DATA --------")
#1.Obtain ID_compound-------------------

#remove non valid inchikeys
notValidInchikeys <- is.na(metadata$inchikey) | nchar(metadata$inchikey) < 10
metadata[notValidInchikeys, "inchikey"] <- NA_character_

if(removeRedCompounds){
  message("\nRedundant metabolites------")
  # IF INCHIKEY!=NA every unique inchikey means a different ID_compound
  #ELSE for every metabolite a different ID_compound

  comparableInchikeys <- nchar(metadata$inchikey) == 27 & !is.na(metadata$inchikey)
  # "IT HAS inchikey" CASE
  # inchies_fact <- as.factor(metadata$inchikey)
  # metadata$ID_compound[!is.na(metadata$inchikey)] <- as.numeric(inchies_fact[!is.na(metadata$inchikey)])
  inchies_fact <- as.factor(metadata$inchikey[comparableInchikeys])
  metadata$ID_compound[comparableInchikeys] <- as.numeric(inchies_fact)

  #"IT HAS NOT inchikey" CASE
  metadata$ID_compound[!comparableInchikeys] <- max(metadata$ID_compound, na.rm = TRUE) + seq_len(sum(!comparableInchikeys))
  metadata$ID_compound <- as.integer(metadata$ID_compound)
  message(glue::glue("
{numCmp - length(unique(metadata$ID_compound))} out of {numCmp} \\
metabolites are redundant and will be eliminated
                   "))
}else{
  metadata$ID_compound <- as.integer(seq_len(nrow(metadata)))
}



#2. Identify repeated spectra----------------
message("\nRedundant spectra------")
metadata$ID_spectra <- NA
#Look for redundant spectra
list_posEspectresduplicats <- .redundantSpectra(listfrag = fragm)
if(length(list_posEspectresduplicats) == 0){
  initVal <- 0
}else{
  #NOMES ELIMINAREM ESPECTRES IGUALS SI TENEN var_espectrals_distintives IDENTIQUES.
  #Pq no volem que quan filtrem per una d'aquestes variables, no hi sigui un espectre que hem eliminat per duplicitat
  #La funcio considera_varespdist() funciona amb ID_spectras i aqui encara no els tenim. Fem servir les posicions com a pseudo ID_spectras
  #substituim posicions pels seus ID_spectras (pq encara no els hem posat)
  message("Analyzing distinctive spectral variables")

  list_posEspectresduplicats_VED <- .applyDistintVars(
    listidespectre_EspDupli = list_posEspectresduplicats,
    dfmetaespect_VarEspDist =
      data.frame(metadata[colnames(metadata) %in% DISTSPECTRALVARS],
                 ID_spectra = seq_len(nrow(metadata))
    ))

  #Assign same ID_spectra to repeated spectra
  for(i in seq_along(list_posEspectresduplicats_VED)){
    metadata$ID_spectra[list_posEspectresduplicats_VED[[i]]] <- i
  }
  initVal <- max(metadata$ID_spectra, na.rm = TRUE)
}

# Assign different ID_spectra to spectra not duplicated (beginning with the last ID_spectra assigned on the last chunk of code)
toassign <- is.na(metadata$ID_spectra)
metadata$ID_spectra[toassign] <- initVal + seq_along(which(toassign))
metadata$ID_spectra <- as.integer(metadata$ID_spectra)
fragments <- list(ID_spectra = metadata$ID_spectra, spectra = fragm)
rm(fragm)
message(glue::glue("
{numMs - length(unique(metadata$ID_spectra))} out of {numMs} \\
spectra are redundant and will be eliminated
                   "))
#3. Restructure data into DF-------
metadata$ID_db <- nameDB
#_3.1 DF Relacional Espectre-Metabolit -------------------------
spectraCompounds <- metadata[, c("ID_spectra","ID_compound")]
#Remove duplicated relations Spectra-Metabolite
spectraCompounds <- distinct(spectraCompounds)

#_3.2_DF MetaMetabolits-------------------------
compounds <- metadata[,
                    c("ID_compound", "ID_db",
                      rownames(varsToParse)[varsToParse$type=="metaboliteVar"])]
compounds <- distinct(compounds, ID_compound, .keep_all = TRUE)

#_3.3 DF MetaEspectres-------------------------
spectra <- metadata[ , c("ID_spectra","ID_db",
                         rownames(varsToParse)[varsToParse$type=="spectraVar"])]
spectra <- distinct(spectra, ID_spectra, .keep_all = TRUE)
notRepeatedIDspect <- !duplicated(fragments$ID_spectra)
fragments$ID_spectra <- fragments$ID_spectra[notRepeatedIDspect]
fragments$spectra <- fragments$spectra[notRepeatedIDspect]


#_3.4_DF BD-------------------------
lastupdate <- format(Sys.time(), "%Y%m%d_%H%M%S")
originalDB <- data.frame(ID_db = unique(metadata$ID_db),
                         lastModification = Sys.Date(),
                         lastRawDataUpdate = lastupdate)

return(list(spectraCompounds = spectraCompounds, spectra = spectra,
            compounds = compounds, originalDB = originalDB,
            fragments = fragments))
}

#funcio per veure data.frames quan llençem un error amb stop()
.print_and_capture <- function(x)
{
  paste(capture.output(print(x)), collapse = "\n")
}


#funcio que agrupa els espectres de frag duplicats. Torna una llista, cada element son posicions d'espectres iguals
#listfrag es NOMES una llista de matrius despectres de fragmentacio.

.redundantSpectra <- function(listfrag){
  message("Searching redundant spectra")
#2 steps. First One-fragment spectra. Then the polyfragmented
  #1. Monofragment spectra (18% of HMDB !!)
  numCols <- unlist(sapply(listfrag, function(x) ncol(x)))
  if(any(numCols < 1))
    stop("There are spectra with less than one fragment")
  #spectra with only one frag
  listfrag_n1 <- lapply(list(seq_along(listfrag), listfrag),
                            function(x) x[numCols == 1])
  mcharge <- sapply(listfrag_n1[[2]], function(x) x["mass-charge", 1])
  mchargeRedund <- unique(mcharge[duplicated(mcharge)])

  #group spectra with the same masscharge (because that means they're equals)
  listfrag_n1 <- lapply(mchargeRedund,
                        function(x) listfrag_n1[[1]][which(mcharge == x)])

  #2. Polyfragmented spectra
  listfrag_n2 <- lapply(list(seq_along(listfrag), listfrag),
                        function(x) x[numCols > 1])
  #sort columns according their mass
  listfrag_n2[[2]] <- lapply(listfrag_n2[[2]], function(x) x[, order(x[1, ])])

  #2A. Fast aproximation: Obtain a fast hash. Only 2 spectra with the same hash can be candidates to be the same spectrum
  hash <- unlist(lapply(listfrag_n2[[2]], function(x)
    sum(x["mass-charge", ]) + (sum(x["intensity",])/max(x["intensity",]))
    ))
  hashRepeatd <- unique(hash[duplicated(hash)])
  if(length(hashRepeatd) > 0){
  #2.B Group spectra with the same hash, which are candidates to be equal. Then compare spectra of the same group
    sublistIdentFragm <- vector("list", length = 0)
    #PARALELITZAR
    listfrag_n2 <- pbapply::pblapply(hashRepeatd, function(hash_repe) {
      #positions of spectra with same hash
      pos_SpctrSameHash <- which(hash == hash_repe)
      SpctrSameHash <- listfrag_n2[[2]][pos_SpctrSameHash]
      pos2Check <- seq_along(pos_SpctrSameHash)
      leftovers <- rep(TRUE, length(pos2Check))
      # A is always the first of the remaining matrices to be compared (leftovers)
      posMatriu_1 <- 1
      while(sum(leftovers) > 1){
        #positions with common spectrum
        commSpctr <- vector("integer", length = 0)
        #Here we compare, of the remaining matrices to be compared, the first matrix with the rest. If a matrix remains, it is no longer compared because it does not match; therefore, there are single matrices that do not even appear in the final list. There are also single matrices that appear in the final list but ONLY in a vector.
        for(posMatriu_2 in (posMatriu_1 + 1):sum(leftovers)){
          comp <- all.equal(SpctrSameHash[[pos2Check[leftovers][posMatriu_1]]],
                            SpctrSameHash[[pos2Check[leftovers][posMatriu_2]]])
          #If the're equals, we keep the 2nd matrix
          if(isTRUE(comp)) commSpctr <- c(commSpctr, posMatriu_2)
        }
        #We keep the first matrix, which we used to compare
        commSpctr <- c(posMatriu_1, commSpctr)
        #we keep the common spectrum in a list
        identFragm <-
          listfrag_n2[[1]][pos_SpctrSameHash][pos2Check][leftovers][commSpctr]
        sublistIdentFragm <- c(sublistIdentFragm, list(identFragm))
        # mark as already completed the identical spectra
        leftovers[leftovers][commSpctr] <- FALSE
      }
      return(sublistIdentFragm)
    })
    listfrag_n2 <- unlist(listfrag_n2, recursive = FALSE)
    #remove groups wth only one element, because that means they're not redundant. Remember here we have not captured all the unique spectra
    listfrag_n2 <- listfrag_n2[sapply(listfrag_n2, function(x) length(x) > 1)]
  }else{
    listfrag_n2 <- list()
  }

  #MERGE redundant spectra
  listidentFragm <- c(listfrag_n1, listfrag_n2)
  return(listidentFragm)
}


#Donada una llista espdupli on cada item son vectors d'idespectres d'espectres duplicats, obtenim una llista on cada item son idespectres despectres duplicats i que tenen les mateixes variables del df dfmetaespect_VarEspDist. Aquest dfmetaespect_VarEspDist es un df amb una columna dfmetaespect_VarEspDist$idespectre i la resta de columnes son de les Variables espectrals distintives.
.applyDistintVars <- function(listidespectre_EspDupli, dfmetaespect_VarEspDist){
  #per cadascun dels llistats amb les posicions d'espectres iguals
  listidespectre_EspDupli_ambsublists <- pbapply::pblapply(listidespectre_EspDupli, function(y) {
    #busquem qun valor tenen les var_espectrals_distintives de cada espectre repetit
    var_esp_dist <- dfmetaespect_VarEspDist[match(y,dfmetaespect_VarEspDist$idespectre),]
    var_esp_dist$idespectre <- NULL
    #ajuntem totes les var_espectrals_distintives d'un mateix row en una amalgama, treiem espais i fiquem minuscues
    a <- tolower(gsub("[[:space:]]", "", do.call("paste", var_esp_dist))) # note the double square brackets)
    #FEm subllistes amb els que comparteixin amalgama
    lapply(unique(a), function(amalgamaunica) {
      y[a==amalgamaunica]
    })
  })
  #Passem les subllistes a llistes
  listidespectre_EspDupli <- unlist(listidespectre_EspDupli_ambsublists, recursive=FALSE)
  #si han quedat llistes amb un sol espectre, les eliminem (ja no estan repetides)
  listidespectre_EspDupli <- listidespectre_EspDupli[vapply(listidespectre_EspDupli, function(x) length(x)>1, FUN.VALUE =TRUE )]
  return(listidespectre_EspDupli)
}

.removeRedundantCompounds <- function(dfmetabolits, dfEspectreMetabolit){
  nummetabolitsoriginals <- dim(dfmetabolits)[1]
  #busquem els valors inchi duplicats per saber quins son els ID_compound duplicats
  valorsinchi_duplicats <- unique(dfmetabolits$inchikey[duplicated(dfmetabolits$inchikey)])
  #Treiem el cas de que no hi hagi inchi. En cas contrari identificaria com iguals
  valorsinchi_duplicats <- valorsinchi_duplicats[valorsinchi_duplicats!="" & valorsinchi_duplicats!="000000-00-0"]

  #trobem el grups de ID_compound iguals
  list_EqualCompoundsId <- lapply(valorsinchi_duplicats, function(x) dfmetabolits$ID_compound[which(dfmetabolits$inchikey==x)])

  if(length(list_EqualCompoundsId)==0){
    message("NO redundant metabolites found")
    return()
  }

  #dels ID_compound iguals, escollim un de cada grup per salvarse
  idmetab_salvats <- vapply(list_EqualCompoundsId, function(x) x[1], FUN.VALUE = 1)
  #treiem el metabolit que es salvara del llista de metabolits a eliminar
  list_EqualCompoundsId<-lapply(list_EqualCompoundsId, function(x) x[-1])

  #A dfEspectreMetabolits SUBSTITUIM ELS IDMETABOILTS identics pel ID_compound salvat
  #1er mirem quines posicions son
  posicions_acanviar <- lapply(list_EqualCompoundsId, function(x) which(dfEspectreMetabolit$ID_compound %in% x) )
  #2on substituim
  for(i in 1:length(posicions_acanviar)){
    dfEspectreMetabolit$ID_compound[posicions_acanviar[[i]]]<-idmetab_salvats[i]
  }
  #3er esborrem els duplicats
  dfEspectreMetabolit <- distinct(dfEspectreMetabolit)

  #Esborrem els ID_compound duplicats de la df_metametabolit
  dfmetabolits <- dfmetabolits[!(dfmetabolits$ID_compound %in% unlist(list_EqualCompoundsId)),]

  message(paste0(length(unlist(list_EqualCompoundsId))," of ",nummetabolitsoriginals," metabolites have been removed due to they share the same inchikey"))
  return(list(dfmetabolits, dfEspectreMetabolit))
}
