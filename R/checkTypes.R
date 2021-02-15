#' check value types
#'
#' check that values meet the types
#'
#' @param values variable values
#' @param types character vector with variable types. names(types) must contain
#' the variable names.
#'
#' @return
#' @noRd

.checkTypes <- function(values, types){
    b <- match(names(values), names(types))
    for(i in which(!is.na(b))){
        .checkType(values[[i]], types[b[i]])
    }
}

#' check type of one value
#'
#' check that value meets the type
#'
#' @param value character defining the variable name.
#' @param type character defining the variable type. name(type) must be the
#'  name of the variable.
#'
#' @return
#' @noRd

.checkType <- function(value, type){
    if(is.name(value)){ #when argument has no value
        return()
    }else if(type=="integer" & is.numeric(value)){
        if(as.integer(value) != value) pchunk <- "integer"
    }else if(type=="integer"){
        pchunk <- "integer"
    }else if(type=="numeric" & !is.numeric(value)){
        pchunk <- "integer or numeric"
    }else if(type=="logical" & !is.logical(value)){
        pchunk <- "logical (TRUE or FALSE)"
    }else if(type=="character" & !is.character(value)){
        pchunk <- "character"
    }
    if(exists("pchunk")){
        stop(paste0("'", names(type),"' argument is expected to be ", pchunk))
    }
}
