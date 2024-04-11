```R
library('GenomeTelescope')

# data from Nadege
# they are supposidely the same species, but they have very different predicted genome sizes (by 100Mbp); So here we compare the two datasets (for convinience named by the prediction small and large respectively)
spectra_small <- read.table('illumina.histo', col.names = c('cov', 'freq'))
spectra_large <- read.table('ont.histo', col.names = c('cov', 'freq'))
# k 
cov_range <- 70:700
coverage_barplot(spectra_small$freq[cov_range], spectra_large$cov[cov_range])

gs_small <- genomescope(spectra_small, '.', k = 21, readlength = 10000)
gs_large <- genomescope(spectra_large, '.', k = 27, readlength = 10000)

gs_small_coefs <- gs_small[[1]]$m$getEnv()
gs_small_1n_cov <- gs_small_coefs$kmercov

gs_large_coefs <- gs_large[[1]]$m$getEnv()
gs_large_1n_cov <- gs_large_coefs$kmercov

# naive genome size est
sum(as.numeric(spectra_small$freq) * as.numeric(spectra_small$cov)) / (gs_small_1n_cov * 2) / 1e6
# 134.7719 Mbp
sum(as.numeric(spectra_large$freq) * as.numeric(spectra_large$cov)) / (gs_large_1n_cov * 2) / 1e6
# 339.0482

# need to remove left residuals
spectra_small$error_subtracted <- spectra_small$freq
spectra_small$error_subtracted[1:200] <- predict(gs_small[[1]], data.frame(x = c(1:200)))
 
coverage_barplot(spectra_small$error_subtracted[1:1000], spectra_small$cov[1:1000])
sum(as.numeric(spectra_small$error_subtracted) * as.numeric(spectra_small$cov)) / (gs_small_1n_cov * 2) / 1e6
# 119.4275 Mbp

# need to remove left residuals for the spectra_large too
spectra_large$error_subtracted <- spectra_large$freq
spectra_large$error_subtracted[1:40] <- predict(gs_large[[1]], data.frame(x = c(1:40)))
  
coverage_barplot(spectra_large$error_subtracted[1:200], spectra_large$cov[1:200])
sum(as.numeric(spectra_large$error_subtracted) * as.numeric(spectra_large$cov)) / (gs_large_1n_cov * 2) / 1e6
# 223.9726 Mbp

genomic_copy_number_small <- round(as.numeric(spectra_small$cov) / (gs_small_1n_cov ))
sum(genomic_copy_number_small * as.numeric(spectra_small$error_subtracted)) / (2 * 1e6)
# 120.93 -> good enough

genomic_copy_number_large <- round(as.numeric(spectra_large$cov) / gs_large_1n_cov )
sum(genomic_copy_number_large * as.numeric(spectra_large$error_subtracted)) / (2 * 1e6)
# 226.2834 -> that's good enough!

cnft_large <- data.frame(copy_number = genomic_copy_number_large, freq = spectra_large$error_subtracted)
cnft_small <- data.frame(copy_number = genomic_copy_number_small, freq = round(spectra_small$error_subtracted))

small <- aggregate(freq ~ copy_number, data = cnft_small, sum)
large <- aggregate(freq ~ copy_number, data = cnft_large, sum)

sum(large[1:21, 'copy_number'] * large[1:21, 'freq']) / (2 * 1e6)
# 221.0249 -> that's not really explaining the difference (if we look at copy-number only up to 21 occurances)

small_Mpb_per_copynumber <- small[, 1] * small[, 2] / 1e6
large_Mbp_per_copynumber <- large[, 1] * large[, 2] / 1e6

head(small)
coverage_barplot(large[,'freq'],)
font_size <- 1
large[, 'copy_number']
width <- 0.4

plot(large_Mbp_per_copynumber, type = "n", xlab = "Copy number", ylab = "Mbp in a perfect uncollapsed representation", 
    ylim = c(0, max(large_Mbp_per_copynumber)), xlim = c(0,15), 
    cex.lab = font_size, cex.axis = font_size, cex.main = font_size, 
    cex.sub = font_size)
    for (i in 1:15) {
        rect(large[i, 1] - width, 0, large[i, 1] + 
            width, large_Mbp_per_copynumber[i], col = "deepskyblue", border = F)
    }
    for (i in 1:15) {
        rect(small[i, 1] - width, 0, small[i, 1] + 
            width, small_Mpb_per_copynumber[i], col = "red", border = F)
    }
    legend('topright', legend = paste(round(c(sum(large_Mbp_per_copynumber), sum(small_Mpb_per_copynumber)), 1), 'Mbp', c('(ONT)', '(Illumina)')), title = 'total sum', pch = 18,  col = c("deepskyblue", "red"), bty = 'n')
    
```