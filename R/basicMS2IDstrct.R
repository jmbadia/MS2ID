#GOAL
#On every parsing algorithm (no matter from what database), once achieved a common data structure (metadata and fragm argument), apply subsequently a common procedure

#ARGUMENTS
#metadata => dataframe containing spectra & metabolite metadata (variable = columns). No id variables must be included (ID_spectra, ID_compound, ID_db)
#fragm => list where each item is a MS2 spectra. Row order must match with metadata row order
#varsToParse => dataframe with the MS2ID variable names as rownames, the BD original names and the type. Must match with the columns of the metadata
#nameDB => name of the DB we are parsing

.basicMS2IDstrct <- function(metadata, fragm, varsToParse, nameDB,
                              lastRawDataUpdate) {
  #rename variables with MS2ID names
  MS2IDname <- match(colnames(metadata), varsToParse$originalNames)
  colnames(metadata)[!is.na(MS2IDname)] <- varsToParse$MS2IDname[MS2IDname[!is.na(MS2IDname)]]
  if(any(is.na(MS2IDname)))
    warning(paste("Following metadata variables are not considered in MS2ID",
                  "and they will be discarted: ",
                  paste(colnames(metadata)[is.na(MS2IDname)], collapse = ", ")))

  if(!all(varsToParse$MS2IDname %in% varsToParse$originalNames))
    warning(paste("The following MS2ID variables have not been found in",
                "metadata argument and they will not be included: ",
                paste(varsToParse$MS2IDname
                      [!varsToParse$MS2IDname %in% colnames(metadata)],
                      collapse =", ")))

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
  message("\nChecking & repairing data------\n")
#remove entries with NA spectra
NASpectra <- is.na(fragm)
metadata <- metadata[!NASpectra,]
fragm <- fragm[!NASpectra]
#assign NA to empty entries
metadata[metadata==""] <- NA_character_
# Count NA / variable
na_count <- sapply(metadata, function(y) sum(length(which(is.na(y)))))
showNA <- paste(round(100*na_count/nrow(metadata)),"%")# % each variable's NA)
names(showNA) <- names(na_count)
message(paste("NA presence on every variable:\n", .print_and_capture(showNA)))

#1.Obtain ID_compound-------------------
message(paste("\nIdentifying repeated metabolites------\n"))
# IF INCHIKEY!=NA every unique inchikey means a different ID_compound
#ELSE for every metabolite a different ID_compound

#remove non valid inchikeys
metadata[(is.na(metadata$inchikey) | metadata$inchikey=="" | metadata$inchikey=="000000-00-0"),"inchikey"] <- NA_character_

# "IT HAS inchikey" CASE
inchies_fact <- as.factor(metadata$inchikey)
metadata$ID_compound[!is.na(metadata$inchikey)] <- as.numeric(inchies_fact[!is.na(metadata$inchikey)])

#"IT HAS NOT inchikey" CASE
metadata$ID_compound[is.na(metadata$inchikey)] <- max(metadata$ID_compound, na.rm = TRUE)+seq_len(sum(is.na(metadata$inchikey)))

#2. Identify repeated spectra----------------
message("\nIdentifying repeated spectra------\n")

#Look for duplicated spectra
list_posEspectresduplicats <- .redundantSpectra(listfrag = fragm)
#Dels 104427 espectres, 12189 son redundants

#NOMES ELIMINAREM ESPECTRES IGUALS SI TENEN var_espectrals_distintives IDENTIQUES.
#Pq no volem que quan filtrem per una d'aquestes variables, no hi sigui un espectre que hem eliminat per duplicitat
#La funcio considera_varespdist() funciona amb ID_spectras i aqui encara no els tenim. Fem servir les posicions com a pseudo ID_spectras
#substituim posicions pels seus ID_spectras (pq encara no els hem posat)
message("Analyzing distinctive spectral variables\n")

list_posEspectresduplicats_VED <- .applyDistintVars(
  listidespectre_EspDupli = list_posEspectresduplicats,
  dfmetaespect_VarEspDist =
    data.frame(metadata[colnames(metadata) %in% DISTSPECTRALVARS],
               ID_spectra = seq_len(nrow(metadata))
  ))

#Assign same ID_spectra to repeated spectra
metadata$ID_spectra <- NA
for(i in seq_along(list_posEspectresduplicats_VED)){
  metadata$ID_spectra[list_posEspectresduplicats_VED[[i]]] <- i
}

# Assign different ID_spectra to spectra not duplicated (beginning with the last ID_spectra assigned on the last chunk of code)
toassign <- is.na(metadata$ID_spectra)
initVal <- max(metadata$ID_spectra, na.rm = TRUE)
metadata$ID_spectra[toassign] <- initVal + seq_along(which(toassign))
metadata$ID_spectra <- as.integer(metadata$ID_spectra)
fragments <- list(ID_spectra = metadata$ID_spectra, spectra = fragm)
rm(fragm)
metadata$ID_compound <- as.integer(metadata$ID_compound)

#3. Restructure data into DF-------
message("\nRestructuring and sieving data ------\n")
metadata$ID_db <- nameDB
#_3.1 DF Relacional Espectre-Metabolit -------------------------
spectraCompounds <- metadata[, c("ID_spectra","ID_compound")]
temporalvalue <- nrow(spectraCompounds)
#Remove duplicated relations Spectra-Metabolite
spectraCompounds <- distinct(spectraCompounds)
message("Relational dataframe spectraCompounds goes from ", temporalvalue, " rows to ", nrow(spectraCompounds) )

#_3.2_DF MetaMetabolits-------------------------
compounds <- metadata[,
                    c("ID_compound", "ID_db",
                      rownames(varsToParse)[varsToParse$type=="metaboliteVar"])]
num_metabolits_original <- nrow(compounds)
compounds <- distinct(compounds, ID_compound, .keep_all = TRUE)

message(paste("Due to duplicated inchikeys, we reduced metabolites number from ",num_metabolits_original,"to",nrow(compounds),"rows"))

#_3.3 DF MetaEspectres-------------------------
spectra<-metadata[,c("ID_spectra","ID_db",rownames(varsToParse)[varsToParse$type=="spectraVar"])]
num_espectres_original <- nrow(spectra)
spectra <- distinct(spectra, ID_spectra, .keep_all = TRUE)
notRepeatedIDspect <- !duplicated(fragments$ID_spectra)
fragments$ID_spectra <- fragments$ID_spectra[notRepeatedIDspect]
fragments$spectra <- fragments$spectra[notRepeatedIDspect]

message(paste("Due to having identical spectras & var_espectrals_distintives, we reduced the spectra number from",num_espectres_original,"to", nrow(spectra),"rows"))

#_3.4_DF BD-------------------------
lastupdate <- format(Sys.time(), "%Y%m%d_%H%M%S")
originalDB <- data.frame(ID_db= unique(metadata$ID_db), lastModification=Sys.Date(), lastRawDataUpdate = lastupdate)

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
  #ho fem en dos temps. Primer els espectres amb una massa8 es mes facil) i despres els altres
  #
  #_____1.8.1. Identif. MonoEspectres duplicats--------------
  #TRACTEM ELSS FRAGMENTS QUE NOMES TENEN UNA MASSA (ES EL 18% de la HUMAN !!)
  message("Searching repeated monospectra")
  numcolumnes<-unlist(sapply(listfrag, function(x) ncol(x)))

  if(any(numcolumnes<1)) stop("Hi han espectres de fragmentacio amb menys d'una massa ")
  #seleccinem els fragments de nomes una massa
  listfragments_n1<- lapply(list(seq_along(listfrag),listfrag), function(x) x[numcolumnes==1])

  #mirem quins espectres tenen la massa que ja ha sortit a un altra espectre
  valors_masscharge <- sapply(listfragments_n1[[2]], function(x) x["mass-charge",1])
  valorsmassa_duplicats <- unique(valors_masscharge[duplicated(valors_masscharge)])

  #agrupem els fragments que tenen la mateixa massa unica de l'espectre, PQ AIXO VOL DIR Q SON ESPECTRES IDENTICS
  listfragmentsidentics_n1 <- lapply(valorsmassa_duplicats, function(x) listfragments_n1[[1]][which(valors_masscharge==x)])

  #_____1.8.2. Identif. NO MonoEspectres duplicats-------------------
  message("Searching repeated polispectra")
  listfragments_n2<- lapply(list(seq_along(listfrag),listfrag), function(x) x[numcolumnes>1])
  #Si no ho hem fet abans, ordenem les columnes de les matrius de fragments per la massa
  listfragments_n2[[2]] <- lapply(listfragments_n2[[2]], function(x) x[,order(x[1,])])

  #A. Primera APROX: Calculem un hash ràpid. Si dos espectres tenen el mateix hash son candidats a ser el mateix espectre
  hash <- unlist(lapply(listfragments_n2[[2]], function(x) sum(x["mass-charge",]) + (sum(x["intensity",])/max(x["intensity",]))))
  hash_repetits <- unique(hash[duplicated(hash)])

  if(length(hash_repetits)>0){#si hi han hash repetits

    #B Les matrius que tenen un mateix hash fan un grup. Les matrius dun mateix grup (hash) son candidates a ser iguals entre elles. Anem a comparar les matrius dins de cada grup

    sublistfragmentsidentics <- vector("list", length=0)

    #PARALELITZAR
    listfragmentsidentics_n2 <- pbapply::pblapply(
      hash_repetits, function(hash_repe) {
      #posicio de les matrius amb aqst mateix hash
      pos_matriusMateixHash<-which(hash == hash_repe)
      #matrius amb aqst mateix hash
      matriusMateixHash<-listfragments_n2[[2]][pos_matriusMateixHash]
      posicionsarevisar <- 1:length(pos_matriusMateixHash)
      #quines posicionsarevisar queden per comparar
      leftovers <- rep(TRUE,length(posicionsarevisar))
      #La matriu A sempre es la primera de les matrius que queden per comparar (leftovers)
      posMatriu_1 <- 1

      while(sum(leftovers)>1){#mentre quedi mes duna matriu per comparar
        #aqui desarem quines posicions ens han sortit q tenen matrius iguals
        coincidencies <- vector("integer", length=0)
        #aqui comparem, de les matrius del grup q FALTEN per comparar, la primera matriu amb la resta. Si queda una matriu, ja no es compara pq no cal; per tant, hi han matrius uniques que ja ni apareixen al llistat final. Tbe hi han matrius uniques q apareixen al llistat final pero SOLES en un vector

        for(posMatriu_2 in (posMatriu_1+1):sum(leftovers)){
          comp <- all.equal(matriusMateixHash[[posicionsarevisar[leftovers][posMatriu_1]]],matriusMateixHash[[posicionsarevisar[leftovers][posMatriu_2]]])
          #Si son iguals, desem la segona matriu
          if(isTRUE(comp)) coincidencies <- c(coincidencies,posMatriu_2)
        }
        #si hem acabat, desem la primera matriu q hem fet servir en la comparacio
        coincidencies <- c(posMatriu_1,coincidencies)
        #desem els fragments de les matrius trobades iguals
        fragmentsidentics<-listfragments_n2[[1]][pos_matriusMateixHash][posicionsarevisar][leftovers][coincidencies]
        #desem en eun elemnt a banda, els fragments anteriors
        sublistfragmentsidentics<-c(sublistfragmentsidentics,list(fragmentsidentics))
        #eliminem les matrius trobades iguals dels candidats a continuar comparant
        leftovers[leftovers][coincidencies]<-FALSE
      }
      return(sublistfragmentsidentics)
    })
  }
  listfragmentsidentics_n2<-unlist(listfragmentsidentics_n2, recursive=FALSE)
  #eliminem els que fan grups de un, perq vol dir q no son espectres duplicats. Recorda que no tots els unics s'han capturat
  listfragmentsidentics_n2 <- listfragmentsidentics_n2[sapply(listfragmentsidentics_n2, function(x) length(x)>1)]

  #ajunten es llistes de fragments identics
  listfragmentsidentics <- c(listfragmentsidentics_n1,listfragmentsidentics_n2)
  #EL num de espectres que esborrarem es...
  message("Dels ",length(listfrag)," espectres, ",length(unlist(listfragmentsidentics))-length(listfragmentsidentics)," son redundants")
  return(listfragmentsidentics)
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
