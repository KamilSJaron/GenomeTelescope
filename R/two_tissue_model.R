#' Fitting two indipendent peaks to the kmer histogram, this model aims to allow for testing of evenness of spacing (expected 1:2)
#'
#' @param x, y vectors with coverages and theuir respective numbers of kmers (frequencies)
#'
#' @param kmerEst, lengthEst, hetEst parameters that serve as initial estimates of kmer coverage, genome length and proportion of kmers in the 1n peak
#'
#' @export

nlsLM_2peak_two_tissue_model <- function(x, y, kmerEst, lengthEst, hetEst = 0.6){
  nlsLM(y ~ ((het       * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
            ((1 - het)  * dnbinom(x, size = kmercov2  / bias, mu = kmercov2))) * length,
        start = list(kmercov = kmerEst, kmercov2 = (2 * kmerEst), bias = 0.5, length = lengthEst, het = hetEst),
        control = list(minFactor=1e-12, maxiter=40))
}

#' @export
predict_1n_peak_nlsLM_2peak_two_tissue_model <- function(model){
  model_env <- model$m$getEnv()
  model_env$het * dnbinom(model_env$x, size = model_env$kmercov   / model_env$bias, mu = model_env$kmercov) * model_env$length
}

#' @export
predict_2n_peak_nlsLM_2peak_two_tissue_model <- function(model){
  model_env <- model$m$getEnv()
  (1 - model_env$het)  * dnbinom(model_env$x, size = model_env$kmercov2  / model_env$bias, mu = model_env$kmercov2) * model_env$length
}

#' PLotting two tissue model
#'
#' @param model to be plotted
#'
#' @export
plot_two_tissue_model <- function(model){
  x <- model$m$getEnv()$x
  y <- model$m$getEnv()$y
  cov_1n <- coef(model)['kmercov']
  cov_2n <- coef(model)['kmercov2']
  monoploid_col <- 'darkorchid4'
  diploid_col <- 'chocolate'

  coverage_barplot(y, x)

  disomic_prediction <- predict_2n_peak_nlsLM_2peak_unconditional(model)
  monosomic_prediction <- predict_1n_peak_nlsLM_2peak_unconditional(model)

  lines(predict(model, response = T) ~ x, lwd = 3)
  lines(monosomic_prediction ~ x, lwd = 3, col = monoploid_col)
  lines(disomic_prediction ~ x, lwd = 3, col = diploid_col)

  lines(c(cov_1n, cov_1n) - 0.5, c(0, max(y)), lty = 2)
  lines(c(cov_2n, cov_2n) - 0.5, c(0, max(y)), lty = 2)

  legend('topright',
         c('kmer histogram','full model', 'monoploid', 'diploid'),
         col = c('deepskyblue','black', monoploid_col, diploid_col),
         lty = c(NA, 1, 1, 1), lwd = 3,
         pch = c(15, NA, NA, NA), bty = 'n')
}
