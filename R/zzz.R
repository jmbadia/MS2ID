#' @noRd
.onLoad <- function(libname, pkgname) {
    names <- c("mzIndexPTR1", "mzIndexPTR2", "mzIndexPTR3",
                          "mzIndexPTR4", "mzIndexPTR5", "mzIndexPTR6")
    assign("MZINDEXPTR", names, envir = asNamespace(pkgname))
    #type of metric (incrm. or decremental)
    decrMetric <- c("topsoe", "squared_chord")
    assign("DECRMETRIC", decrMetric, envir = asNamespace(pkgname))
    incrMetric <- c("cosine", "fidelity", "metricFunc")
    assign("INCRMETRIC", incrMetric, envir = asNamespace(pkgname))

    #Variables to parse and its name in MS2ID or CompoundDB DB ---
    # ROW ORDER will determine COL ORDER in the final MS2ID DB
    mtbolVar_MS2ID <- c("originalCompoundId", "character",
                        "name", "character",
                        "formula", "character",
                        "exactmass", "numeric",
                        "ExactMass", "numeric",
                        "inchikey", "character",
                        "Kegg", "character",
                        "Chebi", "numeric",
                        "PubchemCID", "numeric",
                        "PubchemSID", "numeric",
                        "cas", "character",
                        "smiles", "character")
    #https://github.com/EuracBiomedicalResearch/CompoundDb/issues/62
    mtbolVar_CompoundDB <- c("-", "name", "formula", "exactmass", "-",
                             "inchikey", "-", "-", "-", "-", "cas", "smiles")
    spctraVar_MS2ID <- c("originalSpectrumId", "character",
                         "splash", "character",
                         "adduct", "character",
                         "precursorMz", "numeric",
                         "precursorMassPath",  "character",
                         "collisionEnergy", "character",
                         "polarity", "integer",
                         "msLevel",  "integer",
                         "predicted", "integer",
                         "ionSource",  "character",
                         "fragmentMode", "character",
                         "instrumentType", "character",
                         "instrument", "character",
                         "resolution",  "character",
                         "originalSpectrumDb", "character"
                         )
    spctraVar_CompoundDB <- c(
        "original_spectrum_id", "splash", "adduct", "precursor_mz", "-",
        "collision_energy", "polarity", "ms_level", "predicted", "-", "-",
        "instrument_type", "instrument", "-", "originalSpectrumDb")
    vars2parse <- data.frame(
        MS2IDname = c(mtbolVar_MS2ID, spctraVar_MS2ID)[c(TRUE, FALSE)],
        CompoundDBname = c(mtbolVar_CompoundDB, spctraVar_CompoundDB),
        originalNames = NA,
        MS2IDdataType = c(mtbolVar_MS2ID, spctraVar_MS2ID)[c(FALSE, TRUE)],
        type = c(
            rep("metaboliteVar", length(mtbolVar_MS2ID)/2),
            rep("spectraVar", length(spctraVar_MS2ID)/2)
            ),
        row.names = c(mtbolVar_MS2ID, spctraVar_MS2ID)[c(TRUE, FALSE)]
        )
    assign("VARS2PARSE", vars2parse, envir = asNamespace(pkgname))
    #Variables espectrals distintives: Eliminem els espectres identics  sempre i quan tinguin les V.E.D. (q definim aqui) iguals. PQ? Exemple:Si un espectre esta repetit amb E=20v i E=30v i eliminem el de 30v quan  cerca espectral per E=30v no sortira
    distSpectralVars <- c("collisionEnergy", "polarity", "predicted",
                              "fragmentMode")
    assign("DISTSPECTRALVARS", distSpectralVars, envir = asNamespace(pkgname))
}
