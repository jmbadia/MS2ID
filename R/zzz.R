#' @noRd
.onLoad <- function(libname, pkgname) {
    names <- c("mzIndexPTR1", "mzIndexPTR2", "mzIndexPTR3",
                          "mzIndexPTR4", "mzIndexPTR5", "mzIndexPTR6")
    assign("MZINDEXPTR", names, envir = asNamespace(pkgname))
}
