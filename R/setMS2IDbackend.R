# CREATE MYSQL + bigmemory files-----------------
#library(DBI)/bigmemory/SQL
#overwrite=NA => if dir already exists, neither overwrite it or throw an error.
.setMS2IDbackend <- function(MS2IDdir, DB, overwrite = FALSE){
    if(missing(MS2IDdir))
        stop("'MS2IDdir' is a mandatory argument")
    if(missing(DB))
        stop("'DB' is a mandatory argument")
    if(dir.exists(MS2IDdir)){
        if(!is.na(overwrite)){
            if(overwrite){
                unlink(MS2IDdir, recursive = TRUE)
                dir.create(MS2IDdir)
            }else{
                stop(glue::glue("
                            'MS2IDdir' directory already exists. Use the \\
                            'overwrite' argument if need it"))
            }
        }
    }else{
        dir.create(MS2IDdir)
    }

    #ORDER COLUMNS--------------------
    orderSpectra <- c("ID_spectra",
                      VARS2PARSE$MS2IDname[VARS2PARSE$type == "spectraVar"])
    orderSpectra <- append(orderSpectra, "ID_db", after = 2)
    DB$spectra <- DB$spectra %>%
        relocate(intersect(orderSpectra, names(.)))

    orderCompounds <- c("ID_compound",
                        VARS2PARSE$MS2IDname[VARS2PARSE$type == "metaboliteVar"]
                        )
    orderCompounds <- append(orderCompounds, "ID_db", after = 2)
    DB$compounds <- DB$compounds %>%
        relocate(intersect(orderCompounds, names(.)))


    ## SECTION 1: Create RSQL database in disk file MS2ID.db keeping peakslit apart with ffData
    #SQL database for metadata
    con <- DBI::dbConnect(RSQLite::SQLite(),
                          dbname = file.path(MS2IDdir, "MS2ID.db"))
    s <- sprintf("create table %s(%s, primary key(%s))", "metaSpectrum",
                 paste(names(DB$spectra), collapse = ", "), "ID_spectra")
    DBI::dbExecute(conn = con, statement = s)
    DBI::dbWriteTable(con, "metaSpectrum", DB$spectra, append = TRUE,
                 row.names = FALSE)
    s <- sprintf("create table %s(%s, primary key(%s))", "metaCompound",
                 paste(names(DB$compounds), collapse = ", "), "ID_compound")
    DBI::dbExecute(conn = con, statement = s)
    DBI::dbWriteTable(con, "metaCompound", DB$compounds, append = TRUE,
                 row.names = FALSE)
    DBI::dbWriteTable(con, "crossRef_SpectrComp", DB$spectraCompounds,
                      overwrite=TRUE)
    DBI::dbWriteTable(con, "dbOriginal", DB$originalDB, overwrite=TRUE)

    ## SECTION 2: save mz index-----------
    # 6 dataframes mzIndex_pointer. mzIndex_pointer[[n]] meets mz with the ref
    # spectra that ha such mz on the n most intense peak (menys l'ultim. n=6)
    # on every row, startPos and numSpectra allows to find the
    # IDReference spectra in the mzIndex_raw[ startPos(1:numSpectra) ].
    lastPos <- 0
    for(i in seq_along(DB$mzIdx$spectrIdx)){
        numSpectra <- vapply(DB$mzIdx$spectrIdx[[i]], length, FUN.VALUE = 3)
        startPos <- c(0, cumsum(numSpectra)[-length(numSpectra)]) + lastPos
        mzIndexPTR <- data.frame(id = DB$mzIdx$lstmzIdx[[i]],
                                 numItems = numSpectra, startPos = startPos)
        namePTR <- paste0("mzIndexPTR", i)
        s <- sprintf("create table %s(%s, primary key(%s))", namePTR,
                     paste(names(mzIndexPTR), collapse = ", "), "id")
        DBI::dbExecute(conn = con, statement = s)
        DBI::dbWriteTable(con, namePTR, mzIndexPTR, append = TRUE,
                          row.names = FALSE)
        lastPos <- tail(startPos, 1) + tail(numSpectra ,1)
    }
    # because I can not save a list on SQL but also It is againts SQL nature read
    # pieces of data (from StartPoin to startPoint+numSpectra instead of
    # select-from-where statements) I also use bigmemory
    # Bigmemory file
    mzIndex_body <- unlist(DB$mzIdx$spectrIdx)
    suppressWarnings(
        bigmemory::as.big.matrix(mzIndex_body,
                                 backingfile = "mzIndex_body.bin",
                                 descriptorfile = "mzIndex_body.desc",
                                 backingpath = MS2IDdir)
    )
    #warning normal "Coercing vector to a single-column matrix."

    ## SECTION 3: save SPECTRA-----------
    numPeaks <- unlist(lapply(DB$fragments$spectra, ncol))
    #Fa una taula on s'indica de quin col a quin col hi ha cada fragment list
    df_spectrPTR <- data.frame(id=DB$fragments$ID_spectra,
                               numItems = numPeaks,
                               startPos = c(0, cumsum(numPeaks)[-length(numPeaks)]))
    s <- sprintf("create table %s(%s, primary key(%s))", "spectraPTR",
                 paste(names(df_spectrPTR), collapse = ", "), "id")
    DBI::dbExecute(conn=con, statement = s)
    DBI::dbWriteTable(con, "spectraPTR", df_spectrPTR, append = TRUE,
                      row.names = FALSE)
    DBI::dbListTables(con)
    DBI::dbDisconnect(con)
    #Bigmemory file
    bigmemory::as.big.matrix(do.call(cbind, DB$fragments$spectra),
                             backingfile = "spectra_body.bin",
                             descriptorfile = "spectra_body.desc",
                             backingpath = MS2IDdir)
    message("MS2ID object succesfully created and saved at {MS2IDdir}\n\n")
    invisible(MS2IDdir)
}
