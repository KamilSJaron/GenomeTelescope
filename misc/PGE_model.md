# PGE model

I did something similar to sex-chromosomes, but with the exception I fit there explicitly PGE model - assuming the springtail/sciaridae/cecidomidae PGE type.

It works quite well on BH3-2, I will paste here the testin I have done as well as the functions (should be moved around afterwards).

Exploration (TODO: this was is not yet fixed)

```{R}
# (getting female estimates frst)

### finally, the perfect-most model - explicitely models PGE in the k-mer spectra
PGE_model <- nls_3peak_two_tissue_PGE(x, y, k, estKmercov, lengthEst, rEst, biasEst, fracEst, max_iterations)
PGE_parameters <- coef(PGE_model)

PGE_parameters['kmercov_m']
PGE_parameters['kmercov_p']

lines(predict(PGE_model, response = T) ~ x, col = 'purple', lty = 2, lwd = 3)

plot_nls_3peak_two_tissue_PGE_model(x, y, est_model, estR){
```
