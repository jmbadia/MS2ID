#' Score considering Neutral Loses
#'
#' Score of two MSn spectra considering the neutral Loss between
#' precursor ions. Every query fragment can match with only one reference
#' fragment, considered in the MS2 spectrum (frag_query_x = frag_ref_y) OR in
#' the NL spectra (PrecMass_query - frag_query_x = PrecMass_ref - frag_ref_y =>
#' frag_query_x = frag_ref_y + PrecMass_query - PrecMass_ref = frag_ref_y  -
#' diffPrecMass(precMass_ref-precMass_query).
#' First a matrix with product of intensities where fragments match is obtained.
#' Rows are query fragments, cols reference fragments (match occurs when mass of
#' frag_query_x equals to frag_ref_y OR frag_ref_y - diffPrecMass). Then, using
#' the hungarian algorithm, best convenient fragments pairs (matchs) are
#' paired
#
#' @param spectr1,spectr2 matrix with two numeric rows (m/z and intensity) for every fragment/column.
#' @param massErrorFrag numeric(1) mass error in ppm. Considered in pairing fragments according its m/z.
#' @param diffPrecMass numeric(1) obtained by substracting the precursor masses (diffPrecMass = precMass_reference - precMass_query).
#'
#' @return numeric(1) with the score.
#'
#' @noRd
.NLCosSim <- function(spectr1, spectr2, diffPrecMass, massErrorFrag){
    matrix1 <- .intensityProducts(spectr1, spectr2,
                                 massErrorFrag = massErrorFrag)
    spectr2[1,] <- spectr2[1,] - diffPrecMass
    matrix2 <- .intensityProducts(spectr1, spectr2,
                                 massErrorFrag = massErrorFrag)
    if(!any(matrix2!=0)){#that case NLcos==Cos. CONSIDER avoid rest of the code
        matrix <- matrix1
    }else{
        # keep the max(intensity product) in case (improbable) a queryfrag matchs
        # with the same reference fragment in both matrices
        m <- matrix(c(as.vector(matrix1), as.vector(matrix2)), nrow = 2, byrow=TRUE)
        matrix <- matrix(apply(m, 2 ,max), ncol = ncol(spectr1))
    }
    if(ncol(spectr1) + ncol(spectr2) == 2){#only one query fragment vs one reference fragment
        names(matrix) <- NULL
        return(ifelse(matrix==0, 0, 1))
    }

    if(nrow(matrix) > ncol(matrix)){#solve_LSAP does not work with nrow > ncol
        matrix <- t(matrix)
    }

    optimal <- as.vector(clue::solve_LSAP(matrix, maximum = TRUE))

    cosSimnum <- sum(matrix[cbind(seq_along(optimal), optimal)])
    cosSimden <- sqrt(sum(spectr1[2,]^2) * sum(spectr2[2,]^2))

    return(cosSimnum/cosSimden)
}

#' Matrix of intensity products
#'
# Considering two spectra matrices (m/z, intensity), we obtain a matrix with the
# product of intensities where the fragments match (considering massErrorFrag)
# and 0 if not
#' @param spectr1,spectr2 matrix with two numeric rows (m/z and intensity) for every fragment/column.
#' @param massErrorFrag numeric(1) mass error in ppm. Considered in pairing fragments according its m/z.
#'
#' @return numeric matrix. rows = spectr1 fragments, cols = spectr2 fragments.
#'
#' @noRd
.intensityProducts <- function(spectr1, spectr2, massErrorFrag){
    sapply(seq_len(ncol(spectr1)), function(idFrag) {
        diff <- abs(spectr1[1, idFrag] - spectr2[1,])
        notMatch <- diff > spectr1[1, idFrag] * massErrorFrag/1e6
        diff[notMatch] <- 0
        diff[!notMatch] <- spectr2[2, !notMatch] * spectr1[2, idFrag]
        return(diff)})
}

