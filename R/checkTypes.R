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
    notPresent <- vapply(values, is.null, FUN.VALUE = TRUE) | is.na(values)
    values <- values[!notPresent]
    types <- types[match(names(values), names(types))]
    for(i in seq_along(values)){
        .checkType(values[[i]], types[i])
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
        if(as.integer(value) != value) pchunk <- "type integer"
    }else if(type=="integer"){
        pchunk <- "type integer"
    }else if(type=="numeric" & !is.numeric(value)){
        pchunk <- "type integer or numeric"
    }else if(type=="logical" & !is.logical(value)){
        pchunk <- "type logical (TRUE or FALSE)"
    }else if(type=="dataframe" & !is.data.frame(value)){
        pchunk <- "type data frame"
    }else if(type=="character" & !is.character(value)){
        pchunk <- "type character"
    }else if(type=="MS2ID" & !is(value, "MS2ID")){
        pchunk <- "class MS2ID"
    }else if(type=="Annot" & !is(value, "Annot")){
        pchunk <- "class Annot"
    }else if(type=="function" & !is.function(value)){
        pchunk <- "type function"
    }
    if(exists("pchunk")){
        stop(paste0("'", names(type),"' argument is expected to be of ",
                    pchunk))
    }
}
