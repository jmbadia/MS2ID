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
}
