#' nls_2peak_male_X0 - model for estimating X chromosome size in a species with X0 sex determination
#'
#' @param x - kmer coverages
#'
#' @param y - kmer frequencies
#'
#' @param k - kmer size (defaut 21)
#'
#' @param estKmercov - starting value for k-mer coverage
#'
#' @param estLength - known genome size (estimated from female)
#'
#' @param estR - known heterozygosity (estimated from female)
#'
#' @param max_iterations - maximum iterations that will be done by nls algorithm
#'
#' @export
nls_2peak_male_X0 <- function(x, y, k = 21, estKmercov, estLength, estR, max_iterations = 40){
    model2 = NULL

    cat("trying nls_2peak_kmer_X0 standard algorithm\n")
    # fixed parameters
    r = estR
    length = estLength

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * length * fraction_diploid +
                          dnbinom(x, size = kmercov / bias, mu = kmercov) * length * (1 - fraction_diploid),
                      start = list(kmercov=estKmercov, bias = 0.5, fraction_diploid=0.9),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

#' @export
predict_disomic_portion_2peak_male_X0 <- function(x, r, k, kmercov, bias, length, fraction_diploid){
    ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
    ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * length * fraction_diploid
}

#' @export
predict_monosomic_portion_2peak_male_X0 <- function(x, r, k, kmercov, bias, length, fraction_diploid){
    dnbinom(x, size = kmercov / bias, mu = kmercov) * length * (1 - fraction_diploid)
}
