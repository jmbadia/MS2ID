#' Create a MS2ID backend
#'
#' \code{createMS2ID} processes the data contained in a
#' \code{\link[CompDb-class]{CompDb}} object \insertCite{CompoundDB}{MS2ID} and
#' creates a MS2ID backend by saving data as SQL and bigMemory files in a
#' directory. Later on, that directory will be load as an \linkS4class{MS2ID}
#' object in order to be used repeatedly as reference library with the
#' \code{\link{annotate}} function.
#'
#' @param name character(1) name of the directory where the files will be saved.
#' @param path character(1) with the path where to create the MS2ID backend.
#' @param cmpdb \code{\link[CompDb-class]{CompDb}}(1) object (with MS/MS
#'   spectra).
#' @param noiseThresh A numeric defining the threshold used in the noise
#'   filtering of the MS/MS spectra. e.g. noiseThresh = 0.01 removes peaks with
#'   an intensity less than 1\% of the base peak.
#' @param calcSplash boolean(1) indicating if SPLASH values (when missing) must
#'   be calculated using the \code{\link[splashR]{getSplash}} function.
#' @param calcMmi boolean(1) indicating if compound monoisotopic mass values
#'   (when missing) must be calculated using the
#'   \code{\link[enviPat]{check_chemform}} function.
#' @param overwrite boolean(1) indicating if the function can overwrite
#'   results.
#' @param removeRedComp boolean(1) indicating if elimiate redundant
#'  compounds (based on the inchikey information).
#' @return \code{createMS2ID} returns a character(1) with the MS2ID backend
#'   location. This value must be used as \code{ms2idFolder} parameter in the
#'   \code{MS2ID} constructor.
#' @author Josep M. Badia \email{josepmaria.badia@@urv.cat}
#' @seealso \linkS4class{MS2ID} class.
#' @rdname MS2ID
#' @export
createMS2ID <- function(name = "MS2ID", path = ".", cmpdb, noiseThresh = 0.01,
                        calcSplash = TRUE, calcMmi = TRUE, overwrite = FALSE,
                        removeRedComp = FALSE)
    {
    if(missing(cmpdb))
        stop("Argument 'cmpdb' is required")
    if(!is(cmpdb, "CompDb"))
        stop("Argument 'cmpdb' must be a CompDb object (CompoundDb package)")
    MS2IDdir <- file.path(path, name)
    if(dir.exists(MS2IDdir)){
        if(overwrite){
            unlink(MS2IDdir, recursive = TRUE)
        }else{
            stop(glue::glue("'MS2IDdir' directory already exists. Use \\
            overwrite argument if need it"))
        }
    }

    #.COLLECT() INFO FUNCTION-----------------
    message("\n-------- 1/n. RE-STRUCTURING THE DATA --------")
    src <- src_compdb(cmpdb)
    ## Get a tbl for all tables
    cmp_tbl <- tbl(src, "ms_compound")
    mtdt_tbl <- tbl(src, "metadata")
    msms_tbl <- tbl(src, "msms_spectrum")
    pks_tbl <- tbl(src, "msms_spectrum_peak")
    #syns_tbl <- tbl(src, "synonym")# don't use it now
    ##convert compound & spectrum metadata to real tbl (no sql connection)
    cmp_tbl <- cmp_tbl %>% collect()
    numCmp <- nrow(cmp_tbl)
    msms_tbl <- msms_tbl %>% collect()
    numMs <- nrow(msms_tbl)
    #Merge compounds and spectra into a table
    mrg <- merge(cmp_tbl, msms_tbl, by="compound_id",
                 all.x = FALSE, all.y = FALSE)
    cmpNoSpectra <- sum(!cmp_tbl$compound_id %in% mrg$compound_id)
    if(cmpNoSpectra > 0){
        message(glue::glue("
        {cmpNoSpectra} out of {numCmp} compounds does not have any spectrum \\
        associated and will be eliminated"))
        numCmp <- numCmp - cmpNoSpectra
    }

    #sync order of metadata rows & spectra rows
    #use spectrum_id to sort properly metadata & spectra and remove it
    message("loading spectral data...")
    pks_tbl <- pks_tbl %>%
        select(!peak_id) %>%
        rename("mass-charge" = "mz") %>%
        collect()
    #convert df into a list of peak spectra
    pks_tbl <- split(pks_tbl, f = pks_tbl$spectrum_id)
    message("parsing spectral data")
    fragm_CompDB <- pbapply::pblapply(pks_tbl, function(spctr){
        t(spctr[, c("mass-charge", "intensity") ])
    })

    fragm_spctrID <- vapply(pks_tbl, function(spctr){
        unique(spctr$spectrum_id)
    }, FUN.VALUE = 1)
    varsToParse <- VARS2PARSE
    varsToParse$originalNames <- varsToParse$CompoundDBname
    #order metadata according spectra order
    mrg <- mrg[match(fragm_spctrID, mrg$spectrum_id),]
    if("originalSpectrumDb" %in% names(mrg)){
        nameDB = mrg$originalSpectrumDb
        mrg$originalSpectrumDb <- NULL
    }else{
        nameDB="unknown"
    }
    DB <- .basicMS2IDstrct(metadata = mrg, fragm = fragm_CompDB,
                      varsToParse = varsToParse, nameDB = nameDB,
                      removeRedCompounds = removeRedComp,
                      numCmp = numCmp, numMs = numMs)

    ###TEMPORARILY. we do not comtemplate DB union (check 10 to implement HERE)
    # ALSO elimina metabolits (amb funcio .removeRedundantCompounds()) i spectra duplicats. Només necesari quan juntem diferents

    message("\n-------- 3/n. CALCULATING MISSING VARIABLES --------")
    if(!calcSplash & !calcMmi)
        message("No missing variables to calculate (step avoided by arguments)")
    #Obtain Splash where is not present
    if(calcSplash){
        message("Obtaining spectra's SPLASH values using splashR")
        if("splashR" %in% installed.packages()[, "Package"]){
            notSplash <- is.na(DB$spectra$splash) | DB$spectra$splash %in% c("")
            if(any(notSplash)){
                fragmNoSplash <- match(DB$spectra$ID_spectra[notSplash],
                                       DB$fragments$ID_spectra)
                splashos <- pbapply::pbsapply(
                    DB$fragments$spectra[fragmNoSplash],
                    function(x) splashR::getSplash(t(x))
                    )
                DB$spectra$splash[notSplash] <- splashos
            }

        }else{
            warning(glue::glue(
            "SPLASH values of spectra could not be calculated. The \\
            required package 'splashR' is not installed in the system."))
        }
    }
    #Obtain compound monoisotopicMW using envipat. Only on rows with formula
    if(calcMmi){
        message("Obtaining compound monoisotopic mass using envipat")
        if("enviPat" %in% installed.packages()[, "Package"]){
            noFormula <- is.na(DB$compounds$formula) |
                DB$compounds$formula == ""
            if(!all(noFormula)){
                utils::data(isotopes, package = "enviPat")
                ept <- enviPat::check_chemform(isotopes,
                                               DB$compounds$formula[!noFormula])
                DB$compounds$exactmass[!noFormula][!ept$warning] <-
                    ept$monoisotopic_mass[!ept$warning]
            }
        }else{
            warning(glue::glue(
            "The molecular mass of the compounds could not be calculated. The \\
            required package 'enviPat' is not installed in the system."))
        }
    }
    #spectra treatment
    DB$fragments$spectra <- lapply(DB$fragments$spectra, function(x) {
        # Sort spectra by mass (necessary in order to apply properly COS SIM)
        x[, order(x[1,], decreasing = F), drop = FALSE]

        # Convert any spectrum vector to matrix. Spectrum with only one
        # fragment can be saved as 1x2 vector instead of 2x1 matrix, inducing
        #  mistakes.
        if(is(x, "numeric") & length(x) == 2)
            return(as.matrix(x))
        else if(!is(x, "matrix"))
            stop("Spectral data type unknown")

        #Filter noise
        x[, x["intensity", ] > noiseThresh * max(x["intensity",]), drop=FALSE]
        })

    #order compounds by neutral mass. Need it to identify() (in order to find faster the REF spectra with the neutral mass we want)
    DB$compounds <- DB$compounds[order(DB$compounds$exactmass,
                                       decreasing = FALSE ), ]

    #order by precursor mass the spectra. Need it to identify() (more precisely, in order to find faster the REF spectra with the precurso mass we want)
    DB$spectra <- DB$spectra[order(DB$spectra$precursorMz,
                                   decreasing = FALSE ), ]

    #create mz index  ------------------------------
    message("\n-------- 4/n. OBTAINING THE MZ INDEX --------")
    DB$mzIdx <- .mzIdx(DB)
#saveRDS(DB, "esborraDB20210819.rds")
    #create MYSQL + bigmemory files-----------------
    message("\n-------- 5/n. SAVING THE MS2ID BACKEND FILES --------")
    .setMS2IDbackend(MS2IDdir, DB, overwrite = NA)
    invisible(MS2IDdir)
}
