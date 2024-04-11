nlsLM_2peak_model <- function(x, y, estKmercov, estLength, estR, k = 31){
	nlsLM(y ~ (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
			   ((1-r)^k)         * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * length,
		  start = list(kmercov = estKmercov, bias = 0.5, length = estLength, r = estR),
		  control = list(minFactor=1e-12, maxiter=40))
}

# plot_heterogametic_model_fixed_r_model <- function(het_model){
#        # variable from the model environment
#        x <- het_model$m$getEnv()$x
#        y <- het_model$m$getEnv()$y
#        k <- het_model$m$getEnv()$k
#        r <- het_model$m$getEnv()$r
# 	   predicted_by_model <- predict(het_model, response = T)

#        barplot <- barplot(y ~ x, col = 'deepskyblue', border = F, xlab = 'Coverage', ylab = 'Frequency', ylim = c(0, max(predicted_by_model) * 1.25))
#        lines(predicted_by_model ~ barplot, lwd = 3)

#        estimates <- coef(het_model)
#        # extracting all the fitted values
#        kmercov <-  estimates['kmercov']
#        bias <-  estimates['bias']
#        est_length <- estimates['length']
# 	   heterogameticSize <- estimates['hetSize']

#        disomic_prediction <- predict_disomic_portion_2peak(x, k, kmercov, bias, est_length, heterogameticSize, r)
#        monosomic_prediction <- predict_monosomic_portion_2peak(x, k, kmercov, bias, est_length, heterogameticSize, r)

#        lines(disomic_prediction ~ barplot, lwd = 3, col = 'darkgoldenrod1')
#        lines(monosomic_prediction ~ barplot, lwd = 3, col = 'red')

#        legend('topright',
#               c('kmer histogram','full model', 'autosomes', 'Y + X chromosomes'),
#               col = c('deepskyblue','black', 'darkgoldenrod1', 'red'),
#               lty = c(NA, 1, 1, 1), lwd = 3,
#               pch = c(15, NA, NA, NA), bty = 'n')

#        total_genome <- round(est_length / 1e6, 2)
#        X_chrom_size_est <- round(heterogameticSize / 1e6, 2)
#        het <- round(r * 100, 2)
#        title(paste0("Heterozygosity: ", het , '% X + Y: ', X_chrom_size_est, ' Mbp out of ', total_genome, ' Mbp'))
# }

nlsLM_2peak_model <- function(x, y, estKmercov, estLength, estR, k = 31){
	nlsLM(y ~ (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
			   ((1-r)^k)         * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * length,
		  start = list(kmercov = estKmercov, bias = 0.5, length = estLength, r = estR),
		  control = list(minFactor=1e-12, maxiter=40))
}

predict_2peak_model <- function(model, k = 31){
  model_env <- model$m$getEnv()
    (((2*(1-(1-model_env$r)^k)) * dnbinom(model_env$x, size = model_env$kmercov   / model_env$bias, mu = model_env$kmercov)) +
     ((1-model_env$r)^k)         * dnbinom(model_env$x, size = model_env$kmercov * 2 / model_env$bias, mu = 2 * model_env$kmercov)) * (model_env$length)
}

