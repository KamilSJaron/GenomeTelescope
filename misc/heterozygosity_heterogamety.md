
The idea of heterogametic sex size estimate is that we can estimate the X/Y sizes assuming the heterozygosity from a female sample (or visa reverse for ZW). Here the idea is a bit different. Can we fit 

```R
library('GenomeTelescope')

hoverfly <- read.table('data/idSepSphe2/idSepSphe2.k31.hist.txt', col.names = c('cov', 'freq'))

cov_range <- 4:80
coverage_barplot(hoverfly$freq[cov_range], hoverfly$cov[cov_range])

x <- hoverfly$cov
y <- hoverfly$freq
k <- 31
estKmercov <- 20
estLength <- 400e6
estR <- 0.02
max_iterations <- 20

# basic model of sorts - genomescope style
nls(y ~ ((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
        ((1-r)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2)) * length,
    start = list(r=0, kmercov=estKmercov, bias = 0.5, length=estLength),
    control = list(minFactor=1e-12, maxiter=20))
#         r   kmercov      bias    length 
# 2.138e-02 2.043e+01 3.094e-01 3.877e+08 

genomescope(hoverfly, foldername = 'data/idSepSphe2/genomescope', k = 31, readlength = 15000)
# Model converged het:0.0196 kcov:20.9 err:0.00239 model fit:0.201 len:499780073
#         d         r   kmercov      bias    length 
# 6.123e-02 1.962e-02 2.090e+01 2.011e-01 3.755e+08 
```

Ok, this was the soft-core part. Now we are interested in separating the core genome length from sex-linked genome length. 

From curation: Y is ~72Mb and the suspected X is currently ~69Mb.

```R
curatedHeterogameticSize <- (72 + 69) * 1e6
heterogameticSize <- curatedHeterogameticSize

nlsLM_2peak_heterogametic_model_fixed_size <- function(x, y, estKmercov, estLength, heterogameticSize, estR){
  nlsLM(y ~ (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
            ((1-r)^k)         * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * (length - heterogameticSize) +
            dnbinom(x, size = kmercov   / bias, mu = kmercov) * heterogameticSize,
        start = list(kmercov = estKmercov, bias = 0.5, length = estLength, r = estR),
        control = list(minFactor=1e-12, maxiter=40))
}

het_model <- nlsLM_2peak_heterogametic_model_fixed_size(x[cov_range], y[cov_range], estKmercov, estLength, heterogameticSize, estR)

predict_disomic_portion_2peak <- function(x, k, kmercov, bias, est_length, heterogameticSize, r){
    (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
    ((1-r)^k)         * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * (est_length - heterogameticSize)
}

predict_monosomic_portion_2peak <- function(x, k, kmercov, bias, est_length, heterogameticSize, r){
    dnbinom(x, size = kmercov   / bias, mu = kmercov) * heterogameticSize
}

plot_heterogametic_model_fixed_size_model <- function(het_model){
       # variable from the model environment
       x <- het_model$m$getEnv()$x
       y <- het_model$m$getEnv()$y
       k <- het_model$m$getEnv()$k
       heterogameticSize <- het_model$m$getEnv()$heterogameticSize

       barplot <- barplot(y ~ x, col = 'deepskyblue', border = F, xlab = 'Coverage', ylab = 'Frequency')
       lines(predict(het_model, response = T) ~ barplot, lwd = 3)

       estimates <- coef(het_model)
       # extracting all the fitted values
       kmercov <-  estimates['kmercov']
       bias <-  estimates['bias']
       est_length <- estimates['length']
       r <- estimates['r']

       disomic_prediction <- predict_disomic_portion_2peak(x, k, kmercov, bias, est_length, heterogameticSize, r)
       monosomic_prediction <- predict_monosomic_portion_2peak(x, k, kmercov, bias, est_length, heterogameticSize, r)

       lines(disomic_prediction ~ barplot, lwd = 3, col = 'darkgoldenrod1')
       lines(monosomic_prediction ~ barplot, lwd = 3, col = 'red')

       legend('topright',
              c('kmer histogram','full model', 'autosomes', 'Y + X chromosomes'),
              col = c('deepskyblue','black', 'darkgoldenrod1', 'red'),
              lty = c(NA, 1, 1, 1), lwd = 3,
              pch = c(15, NA, NA, NA), bty = 'n')

       total_genome <- round(est_length / 1e6, 2)
       X_chrom_size_est <- round(heterogameticSize / 1e6, 2)
       het <- round(r * 100, 2)
       title(paste0("Heterozygosity: ", het , '% assumed X+Y size: ', X_chrom_size_est, ' Mbp out of total ', total_genome, ' Mbp'))
}

plot_heterogametic_model_fixed_size_model(het_model)
```

Ok, let's now do the het ~ het estimate.


```R
hetSizes <- c(50e6,100e6, curatedHeterogameticSize, 200e6, 250e6, 300e6)

hoverfly_wrapper <- function(one_hetSize){
    nlsLM_2peak_heterogametic_model_fixed_size(x[cov_range], y[cov_range], estKmercov, estLength, one_hetSize, estR)
}

het_models <- lapply(hetSizes, hoverfly_wrapper)
heterozygosities <- sapply(het_models, function(x){coef(x)['r']})

par(mfrow = c(2, 3))
lapply(het_models, plot_heterogametic_model_fixed_size_model)
dev.off()

plot(heterozygosities * 100, hetSizes / 1e6, ylim = c(0, max(hetSizes) / 1e6), type = 'l', ylab = "Assumed size of X + Y [Mbp]", xlab = 'Estimated heterozygosity [%]')
points((heterozygosities * 100)[3], (hetSizes / 1e6)[3], pch = 20, col = 'purple')
```

Great, all that works, can we do it other way around? Fitting the size while assuming heterozygosity?

```R
nlsLM_2peak_heterogametic_model_fixed_r <- function(x, y, estKmercov, estLength, estHeterogameticSize, r){
  nlsLM(y ~ ((2*(1-(1-r)^k))  * dnbinom(x, size = kmercov   / bias, mu = kmercov) +
             ((1-r)^k)        * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * (length - heterogameticSize) +
                                dnbinom(x, size = kmercov   / bias, mu = kmercov) * heterogameticSize,
        start = list(kmercov = estKmercov, bias = 0.5, length = estLength, heterogameticSize = estHeterogameticSize),
        control = list(minFactor=1e-12, maxiter=40))
}

het_model <- nlsLM_2peak_heterogametic_model_fixed_r(x[cov_range], y[cov_range], estKmercov, estLength, heterogameticSize, estR)

predict_disomic_portion_2peak_fixed_r <- function(x, k, kmercov, bias, est_length, heterogameticSize, r){
    (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
    ((1-r)^k)         * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * (est_length - heterogameticSize)
}

predict_monomic_portion_2peak_fixed_r <- function(x, k, kmercov, bias, est_length, heterogameticSize, r){
    dnbinom(x, size = kmercov   / bias, mu = kmercov) * heterogameticSize
}
```