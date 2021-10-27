#' ESI-MS adducts
#'
#' List of common adducts observed for ESI-MS measurements in soft
#'  positive and negative ionization modes. Obtained from the Envipat package
#'
#' @format A data frame with 49 rows and 9 variables:
#' \describe{
#'   \item{Name}{Adduct name}
#'   \item{calc}{Equation for calculating adduct m/z from uncharged non-adduct molecular mass M (m/z =M/z + X)}
#'   \item{Charge}{z}
#'   \item{Mult}{1/z}
#'   \item{Mass}{X}
#'   \item{Ion_mode}{Ionization mode (1=positive, 0=negative)}
#'   \item{Formula_add}{Adduct chemical formula to be added}
#'   \item{Formula_ded}{Adduct chemical formula to be subtracted}
#'   \item{Multi}{Factor to multiply chemical formula with}
#' }
#' @source \url{https://cran.r-project.org/package=enviPat}
#' @noRd
"adducts"
