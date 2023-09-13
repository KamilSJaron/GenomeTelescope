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
nls_2peak_male_XY <- function(x, y, k, estKmercov, estLength, estR, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_kmer_XY standard algorithm\n")

    # fixed parameter
    r = estR

    try(model2 <- nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
                          ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * disomic_length +
                          dnbinom(x, size = kmercov / bias, mu = kmercov) * monosomic_length,
                      start = list(kmercov=estKmercov, bias = 0.5, monosomic_length = 0.1 * estLength, disomic_length = 0.9 * estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model2)
}

#' @export
predict_disomic_portion_2peak_male_XY <- function(x, r, k, kmercov, bias, disomic_length){
    ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
    ((1-r)^k)        * dnbinom(x, size = kmercov * 2 / bias, mu = kmercov * 2)) * disomic_length
}

#' @export
predict_monosomic_portion_2peak_male_XY <- function(x, kmercov, bias, monosomic_length){
    dnbinom(x, size = kmercov / bias, mu = kmercov) * monosomic_length
}

#' @export
nls_4peak_male_XY <- function(x, y, k, estKmercov, estLength, estR, estD, max_iterations){
    model4 = NULL

    cat("trying nls_4peak_kmer_XY standard algorithm\n")

    # fixed parameter
    r = estR
    d = estD

    try(model4 <- nls(y ~ (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)       +
                          (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)   +
                          (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3)   +
                          (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4))  * disomic_length +
                          ((1 - d)                                                                     * dnbinom(x, size = kmercov   / bias, mu = kmercov)       +
                          (d                                                                           * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2))) * monosomic_length,
                      start = list(kmercov=estKmercov, bias = 0.5, monosomic_length = 0.5 * estLength, disomic_length = 0.5 * estLength),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = TRUE)

    return(model4)
}

#' @export
predict_disomic_portion_4peak_male_XY <- function(x, r, d, k, kmercov, bias, disomic_length){
    (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)       +
    (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)   +
    (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3)   +
    (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4)) * disomic_length
}

#' @export
predict_monosomic_portion_4peak_male_XY <- function(x, d, kmercov, bias, monosomic_length){
    ((1 - d)                                                                     * dnbinom(x, size = kmercov   / bias, mu = kmercov)       +
    (d                                                                           * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2))) * monosomic_length
}

#' @export
plot_XY_model <- function(XY_model){
       # variable from the model environment
       x <- XY_model$m$getEnv()$x
       y <- XY_model$m$getEnv()$y
       k <- XY_model$m$getEnv()$k
       female_r <- XY_model$m$getEnv()$r
       female_d <- XY_model$m$getEnv()$d

       barplot <- barplot(y ~ x, col = 'deepskyblue', border = F, xlab = 'Coverage', ylab = 'Frequency')
       lines(predict(XY_model, response = T) ~ barplot, lwd = 3)

       estimates <- coef(XY_model)
       # extracting all the fitted values
       kmercov <-  estimates['kmercov']
       bias <-  estimates['bias']
       disomic <- estimates['disomic_length']
       monosomic <- estimates['monosomic_length']


       disomic_prediction <- predict_disomic_portion_4peak_male_XY(x, female_r, female_d, k, kmercov, bias, disomic)
       monosomic_prediction <- predict_monosomic_portion_4peak_male_XY(x, female_d, kmercov, bias, monosomic)

       lines(disomic_prediction ~ barplot, lwd = 3, col = 'darkgoldenrod1')
       lines(monosomic_prediction ~ barplot, lwd = 3, col = 'red')

       legend('topright',
              c('kmer histogram','full model', 'autosomes', 'X chromosomes'),
              col = c('deepskyblue','black', 'darkgoldenrod1', 'red'),
              lty = c(NA, 1, 1, 1), lwd = 3,
              pch = c(15, NA, NA, NA), bty = 'n')

       total_genome <- round((monosomic + disomic) / 1e6, 2)
       X_chrom_size_est <- round(monosomic / 1e6, 2)
       title(paste0('Estimated size of X + Y: ', X_chrom_size_est, ' Mbp out of total ', total_genome, ' Mbp'))
}
