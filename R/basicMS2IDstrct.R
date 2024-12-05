#' Adapt the metadata to the MS2ID structure
#'
#' @param metadata dataframe containing spectra & metabolite metadata (variable
#'   = columns). id variables can be included (the name must be compound_id,
#'    spectrum_id and ID_db)
#' @param fragm list where each item is a MS2 spectra. Row ORDER must match with
#'   metadata row ORDER
#' @param varsToParse dataframe with the MS2ID variable names as rownames, the
#'   BD original names and the type. Must match with the columns of the metadata
#' @param nameDB char(1) with the name we want for the DB we are parsing
#' @param numCmp int(1) number of compounds in the original DB
#' @param numMs int(1) number of spectra in the original DB
#' @param removeRedundant character(1) eliminate redundant compounds and/or spectra. Available options: "compounds", "spectra" or "both"
#' @noRd
.basicMS2IDstrct <- function(metadata, fragm, varsToParse,
                             nameDB="unknown", numCmp, numMs,
                             removeRedundant = "compounds")
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
  #we accept character in ID_compound (because it is the compoundb type)
  classOriginal[colnames(metadata) == "ID_compound" &
                    classOriginal == "character"] <- "integer"

  #if it's mandatory numeric (classMS2ID), we accept integer (classOriginal)
  classOriginal[classOriginal == "integer" &
                    classMS2ID == "numeric"] <- "numeric"
  if(!identical(classOriginal, classMS2ID)){
      notEqualClass <- which(classOriginal != classMS2ID)

      messages <- sapply(notEqualClass, function(row) {
          paste0("\n", colnames(metadata)[row], ": ", classOriginal[row],
                  " to ", classMS2ID[row])
      })

      messages <- paste("The following variables do not have the proper data type and must be converted: ",
                        paste(messages, collapse = ", "))
      stop(messages)
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

# compound ID
if ("ID_compound" %in% colnames(metadata)) {# col ID_compound already exists
    # Ensure column A is treated as character (to handle mixed types)
    metadata$ID_compound <- as.character(metadata$ID_compound)
    # Identify non-integer and NA values
    is_integer <- grepl("^[0-9]+$", metadata$ID_compound)

    # must throw an error, otherwhise u're removing info (duplicated "a" would be converted into different id)
    stopifnot("compound id must be integers with no NA's" = sum(is_integer) == nrow(metadata))


    # # Find the maximum integer value in the column
    # max_int <- max(as.numeric(metadata$ID_compound[is_integer]),
    #                na.rm = TRUE)
    # # Replace non-integers with incremental numbers starting from max_int + 1
    # metadata$ID_compound[!is_integer] <- seq(max_int + 1, by = 1,
    #                                          length.out = sum(!is_integer))

    metadata$ID_compound <- as.integer(metadata$ID_compound)
}else{#when no id compound provided, replenish with sequential numbers
    metadata$ID_compound <- as.integer(seq_len(nrow(metadata)))
}

if(removeRedundant %in% c("compounds", "both")){
  message("\nRedundant metabolites------")
  # IF INCHIKEY!=NA every unique inchikey means a different ID_compound
  #ELSE for every metabolite a different ID_compound

  comparableInchikeys <- nchar(metadata$inchikey) == 27 & !is.na(metadata$inchikey)
  # "IT HAS inchikey" CASE
  metadata$ID_compound[comparableInchikeys] <- metadata[comparableInchikeys,] |>
      group_by(inchikey) |>
      mutate(ID_compound = if(n() > 1) min(ID_compound) else ID_compound) |>
      ungroup() |> select(ID_compound) |> unlist()
  # inchies_fact <- as.factor(metadata$inchikey[comparableInchikeys])
  # metadata$ID_compound[comparableInchikeys] <- as.numeric(inchies_fact)

  #"IT HAS NOT inchikey" CASE: dont remove redundant
  # metadata$ID_compound[!comparableInchikeys] <- max(metadata$ID_compound, na.rm = TRUE) + seq_len(sum(!comparableInchikeys))

  metadata$ID_compound <- as.integer(metadata$ID_compound)
  message(glue::glue("
{numCmp - length(unique(metadata$ID_compound))} out of {numCmp} \\
metabolites are redundant and will be eliminated
                   "))
}

#2. Spectra --------
#2.1 Id spectra----------------
# assign id spectra
if ("ID_spectra" %in% colnames(metadata)) {# col ID_spectra already exists
    # Ensure column A is treated as character (to handle mixed types)
    metadata$ID_spectra <- as.character(metadata$ID_spectra)
    # Identify non-integer and NA values
    is_integer <- grepl("^[0-9]+$", metadata$ID_spectra)

    # must throw an error, otherwhise u're removing info (e. g. duplicated "a" would be converted into different id)
    stopifnot("spectral id must be integers with no NA's" = sum(is_integer) == nrow(metadata))

    metadata$ID_spectra <- as.integer(metadata$ID_spectra)
}else{#when no id spectra provided, replenish with sequential numbers
    metadata$ID_spectra <- as.integer(seq_len(nrow(metadata)))
}

#2.2 Identify repeated spectra----------------
if(removeRedundant %in% c("spectra", "both")){
    message("\nRedundant spectra------")
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
      missingVar <- DISTSPECTRALVARS[!DISTSPECTRALVARS %in% colnames(metadata)]
      if(length(missingVar) != 0)
          paste("To remove redundant spectra column/s '%s' must be",
          "present in their metadata in order to differentiate them. Please",
          "add the necessary information or disable the removing spectra",
          "option (i.e. argument removeRedundant to 'compounds' or 'none')") |>
          sprintf(missingVar) |>
          stop()

      list_posEspectresduplicats_VED <- .applyDistintVars(
        listidespectre_EspDupli = list_posEspectresduplicats,
        dfmetaespect_VarEspDist =
          data.frame(metadata[colnames(metadata) %in% DISTSPECTRALVARS],
                     spectra_pos = seq_len(nrow(metadata))
        ))

      #Assign same ID_spectra to repeated spectra
      for(i in seq_along(list_posEspectresduplicats_VED)){
          #among all identical spectra, get position of the first one
          identSpectr_pos <- list_posEspectresduplicats_VED[[i]]
          #use its id with all the obtained identical spectra
          metadata$ID_spectra[identSpectr_pos] <- metadata$ID_spectra[identSpectr_pos[1]]
      }
    }
}

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
colNams <- rownames(varsToParse)[varsToParse$type=="metaboliteVar"] |>
    union(c("ID_compound", "ID_db")) #merge with no duplicities
compounds <- metadata[, colNams]
compounds <- distinct(compounds, ID_compound, .keep_all = TRUE)

#_3.3 DF MetaEspectres-------------------------
colNams <- rownames(varsToParse)[varsToParse$type=="spectraVar"] |>
    union(c("ID_spectra", "ID_db")) #merge with no duplicities
spectra <- metadata[, colNams]

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


#Donada una llista espdupli on cada item son vectors de les posicions d'espectres duplicats, obtenim una llista on cada item son posicions despectres duplicats i que tenen les mateixes variables del df dfmetaespect_VarEspDist. Aquest dfmetaespect_VarEspDist es un df amb una columna dfmetaespect_VarEspDist$spectra_pos i la resta de columnes son de les Variables espectrals distintives.
.applyDistintVars <- function(listidespectre_EspDupli, dfmetaespect_VarEspDist){
    stopifnot("No spectra_pos column found" = "spectra_pos" %in% names(dfmetaespect_VarEspDist),
              "No distinctive variables found" = ncol(dfmetaespect_VarEspDist) > 1)
  #per cadascun dels llistats amb les posicions d'espectres iguals
  listidespectre_EspDupli_ambsublists <- pbapply::pblapply(listidespectre_EspDupli, function(y) {
    #busquem qun valor tenen les var_espectrals_distintives de cada espectre repetit
    var_esp_dist <- dfmetaespect_VarEspDist[match(y, dfmetaespect_VarEspDist$spectra_pos),]
    var_esp_dist$spectra_pos <- NULL
    #ajuntem TOTES les var_espectrals_distintives d'un mateix row en una amalgama (treiem espais i fiquem minuscues)
    a <- tolower(gsub("[[:space:]]", "", do.call("paste", var_esp_dist))) # note the double square brackets)
    #FEm subllistes amb els que comparteixin amalgama
    lapply(unique(a), function(amalgamaunica) {
      y[a == amalgamaunica]
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
