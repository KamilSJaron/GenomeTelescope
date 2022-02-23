#' Given mean +/- stderr, report min and max value within 2 SE
#'
#' @param table, a vector of two element - mean and sd
#'
#' @export

min_max <- function(table){
    ##return (c( abs(table[1]) - 2*abs(table[2]) , abs(table[1])+ 2*abs(table[2])))
    return (c(table[1] - 2*table[2], table[1]+ 2*table[2]))
}

