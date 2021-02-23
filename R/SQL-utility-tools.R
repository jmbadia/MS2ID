## Utility functions to deal with SQL.

#' Append where sentences
#'
#' '.appendSQLwhere' adds a WHERE clause to a vector of WHERE clauses.
#' The aim is to collapse them all - with a 'AND' link - when eventually we
#' construct the SQL sentence
#'
#' @param base character(1) defining the column to apply the condition
#' (e.g. WHERE base BETWEEN condition1 AND condition2)
#' @param condition vector(n) defining conditions. Depends on the mode configuration
#' mode=EQUAL => n=1. E.g. "WHERE base IN=condition"
#' mode=IN => all n conditions collapse. E.g. "WHERE base IN = ('cond1', 'cond2',
#'  'cond3')
#' mode=BETWEEN => n=2. E.g. "WHERE base BETWEEN condition1 AND condition2"
#' @param mode character describing the type of WHERE sentence (EQUAL, IN, BETWEEN )
#' @param whereVector character(n) vector with previous n WHERE clauses
#'
#' @return A vector of characters. Each position is a WHERE sentence
#' @noRd

.appendSQLwhere <- function(base, condition, mode="EQUAL", whereVector){
    if(mode == "EQUAL"){
        if(is.character(condition)) condition <- paste0("'", condition, "'")
        rsltString <- paste(base, "=", condition)
    } else if(mode == "IN"){
        condition <- paste(condition, collapse = ", ")
        rsltString <- paste(base, "IN (", condition,")")
    } else if(mode == "BETWEEN"){
        if(length(condition)!=2) stop("2 conditions are required")
        rsltString <- paste(base, "BETWEEN", condition[1], "AND",
                            condition[2])
    } else stop("The 'mode' choosed is not contemplated")

    if(missing(whereVector)) whereVector <- vector()
    whereVector[length(whereVector)+1] <- rsltString
    return(whereVector)
}

#' Obtain records from a SQL query
#'
#' @param MS2ID MS2ID object containing the SQL database
#' @param select character(n) containing the names of the columns
#'  to read data from
#' @param from character(n) with the names of the tables to be consulted
#' @param where character(n) with n WHERE clauses
#'
#' @return Returns the result of the query as a data frame
#' @noRd
.getSQLrecords <- function(MS2ID, select, from, where){
    where <- ifelse(missing(where),
                    "", paste("WHERE", paste(where, collapse = " AND ")))
    DBI::dbGetQuery(MS2ID@dbcon,
                    paste("SELECT", select, "FROM", from, where))
}
