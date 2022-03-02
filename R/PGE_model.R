#' Fitting two indipendent peaks to the kmer histogram, this model aims to allow for testing of evenness of spacing (expected 1:2)
#'
#' @param x, y vectors with coverages and theuir respective numbers of kmers (frequencies)
#'
#' @param estR - heterozygosity, a fixed parameter treated as a constant
#'
#' @param estKmercov, estLength, biasEst parameters that serve as initial estimates of kmer coverage, genome length and overdispersal
#'
#' @export

nls_3peak_two_tissue_PGE <- function(x, y, k, estKmercov, estLength, estR, biasEst, estFrac, max_iterations){
    model2 = NULL

    cat("trying nls_2peak_two_tissue_constant_shift standard algorithm\n")
    # fixed parameters
    r = estR                 # autosomal heterozygosity
    # length = estLength       # haploid genome length
    # fraction_diploid = estFrac # estimated fraction of Autosomes

    try(model2 <- nls(y ~ ((1-(1-r)^k) * dnbinom(x, size = kmercov_p / bias, mu = kmercov_p) +                                                             # paternal heterozygous autosomes
                           (1-(1-r)^k) * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) +                                                             # maternal heterozygous autosomes
                           ((1-r)^k)   * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid +      # homozygous autosomes
                           1           * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m)                            * length * (1 - fraction_diploid), # maternally inherited X autosomes
                      start = list(kmercov_p=estKmercov, kmercov_m=estKmercov, bias=biasEst, length=estLength, fraction_diploid = estFrac),
                      control = list(minFactor=1e-12, maxiter=max_iterations)), silent = F)

    return(model2)
}

#' @export
predict_maternal_3peak_PGE <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) +
  ((1-r)^k) * 1/2  * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid +
  1                * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m)                            * length * (1 - fraction_diploid)
}

#' @export
predict_maternal_3peak_PGE_het <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m)) * length * fraction_diploid
}

#' @export
predict_maternal_3peak_PGE_autosomes <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) +
  ((1-r)^k) * 1/2  * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid
}

#' @export
predict_maternal_3peak_PGE_X <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
   dnbinom(x, size = kmercov_m / bias, mu = kmercov_m) * length * (1 - fraction_diploid)
}

#' @export
predict_paternal_3peak_PGE <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_p / bias, mu = kmercov_p) +
  ((1-r)^k) * 1/2 * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid
}

#' @export
predict_paternal_3peak_PGE_het <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  ((1-(1-r)^k)     * dnbinom(x, size = kmercov_p / bias, mu = kmercov_p)) * length * fraction_diploid
}

#' @export
predict_shared_3peak_PGE <- function(x, r, k, kmercov_m, kmercov_p, bias, length, fraction_diploid){
  (((1-r)^k) * dnbinom(x, size = (kmercov_p + kmercov_m) / bias, mu = kmercov_p + kmercov_m)) * length * fraction_diploid
}

#' @export
plot_nls_3peak_two_tissue_PGE_model <- function(x, y, est_model, estR){

  fitted_parameters <- coef(est_model)
  cov_m <- fitted_parameters['kmercov_m']
  cov_p <- fitted_parameters['kmercov_p']
  cov <- cov_m + cov_p

  matermal_prediction <- predict_maternal_3peak_PGE(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  matermal_prediction_het <- predict_maternal_3peak_PGE_het(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  matermal_prediction_A <- predict_maternal_3peak_PGE_autosomes(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  matermal_prediction_X <- predict_maternal_3peak_PGE_X(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  paternal_prediction <- predict_paternal_3peak_PGE(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  paternal_prediction_het <- predict_paternal_3peak_PGE_het(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])
  shared_prediction <- predict_shared_3peak_PGE(x, estR, k, cov_m, cov_p, fitted_parameters['bias'], fitted_parameters['length'], fitted_parameters['fraction_diploid'])

  barplot <- barplot(y ~ x, col = 'deepskyblue', border = F, xlab = 'Coverage', ylab = 'Frequency', ylim = c(0, max(y) * 1.1))

  lines(predict(est_model, response = T) ~ barplot, lwd = 3)

  # lines(matermal_prediction ~ barplot, lwd = 3, col = 'darkorange2')
  lines(matermal_prediction_het ~ barplot, lwd = 3, col = 'firebrick3')
  lines(paternal_prediction_het ~ barplot, lwd = 3, col = 'darkblue')
  lines(shared_prediction ~ barplot, lwd = 3, col = 'darkorchid3')
  lines(matermal_prediction_X ~ barplot, lwd = 3, col = "gold")


  # lines(c(cov_m, cov_m), c(-1e10, 1e10), lwd = 2, lty = 2, col = 'firebrick3')
  # lines(c(cov_p, cov_p), c(-1e10, 1e10), lwd = 2, lty = 2, col = 'darkblue')
  # lines(c(cov, cov), c(-1e10, 1e10), lwd = 2, lty = 2, col = 'darkorchid3')

  legend('topright',
         c('kmer histogram','full model', 'shared A', 'maternal X chromosomes', 'maternal heterozygous A', 'paternal heterozygous A'),
         col = c('deepskyblue','black', 'darkorchid3', 'gold', 'firebrick3', 'darkblue'),
         lty = c(NA, 1, 1, 1, 1, 1), lwd = 3,
         pch = c(15, NA, NA, NA, NA, NA), bty = 'n')

  title(paste0('Explicitly fitted PGE model to A. fusca'))

}
