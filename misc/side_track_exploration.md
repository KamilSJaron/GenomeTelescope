## Side track exploration

This is more a dump file, so I can leave somewhere code that I don't want to get completelly lost, but I belibe leads nowhere.

### Totally alternative ways to get 1n and 2n coverage for two tissue model

and some exploration of **kernel smoothing** as replacement.

```{R}
input_file1 <- 'data/springtails/Afus1_k21_truncated.hist' # high coverage male

kmer_spectrum <- read.table(input_file1, col.names = c('coverage', 'frequency'))

range <- 15:150
adjust = 1

second_deriv <- diff(sign(diff(kmer_spectrum$frequency[range])))

peak_covs <- kmer_spectrum$coverage[which(second_deriv == -2) + 1]
peak_heights <- kmer_spectrum$frequency[which(second_deriv == -2) + 1]
# head(peak_covs[order(peak_heights, decreasing=T)])

ks <- density(kmer_spectrum[range, 'coverage'], bw = "SJ", adjust = adjust, weights = kmer_spectrum[range, 'frequency'] / sum(kmer_spectrum[range, 'frequency']))
```

However, it really struggles with cases then the two peaks are largely overlaping. I think non-linear regression model fitting is a lot better way to go


### Basics of Non-linear least sequares

This really early fooling around

```{R}
source('misc/unsorted_genome_models.R')

input_file4 <- 'data/springtails/WW5-3_k21_truncated.hist' # low coverage female
kmer_spectrum <- read.table(input_file4, col.names = c('coverage', 'frequency'))
range <- 4:60

plot(frequency ~ coverage, data = kmer_spectrum[range, ])

x <- kmer_spectrum[range, 'coverage']
y <- kmer_spectrum[range, 'frequency']
estKmercov <- 15
max_iterations <- 40
VERBOSE=1
### original - works
genomescope <- nls_4peak(x, y, 21, estKmercov, sum(y * x) / estKmercov, max_iterations)

### a bit simplified - WORKS
genomescope_2peak <- nls_2peak_kmer_explicit(x, y, 21, estKmercov, sum(y * x) / estKmercov, max_iterations)
lines(predict(genomescope_2peak, response = T) ~ kmer_spectrum$coverage[range])

biasEst <- coef(genomescope_2peak)['bias']
lengthEst <- coef(genomescope_2peak)['length']
kmerCoverageEst<- coef(genomescope_2peak)['kmercov']
rEst <- coef(genomescope_2peak)['r']

conditional_model <- nls_2peak_conditional(x, y, 21, kmerCoverageEst, lengthEst, rEst, biasEst, max_iterations)
lines(predict(conditional_model, response = T) ~ kmer_spectrum$coverage[range], col = 'green')

### simplified - DOES NOT
genome_model <- nls_2peak(x, y, estKmercov, max_iterations)
lines(predict(genome_model, response = T) ~ kmer_spectrum$coverage[range])

# confint2(genome_model, level = 0.95, method = c("asymptotic"))
# est_1n <- coef(genome_model)[1]
# est_2n <- coef(genome_model)[2]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # 1n / 2n kmer coverage estimates to genome coverage
# k = 21
# readlength = 150
# kmer2genome_coverage <- function(est, k, readlength){
# 	return((est * readlength) / (readlength - k + 1))
# }
#
# kmer2genome_coverage(est_1n, k, readlength)
# kmer2genome_coverage(est_2n, k, readlength)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

########
# Quest to find reasonable starting values
########

## First we see what happens when the max peak is the kmercoverage (typically the homozygous peak) for the plot
numofReads = sum(as.numeric(x*y))/(readlength-k+1)
estKmercov  = x[which(y==max(y))][1]
estCoverage1 = estKmercov1*readlength/(readlength-k)
estLength   = numofReads*readlength/estCoverage1
max_iterations = 4
VERBOSE = T
nls1    = nls_4peak(x, y, k, estKmercov1, estLength1, max_iterations)

## Second we half the max kmercoverage (typically the heterozygous peak)
estKmercov2  = estKmercov1 / 2 ##2.5
estCoverage2 = estKmercov2*readlength/(readlength-k)
estLength2   = numofReads*readlength/estCoverage2

nls2 = nls_4peak(x, y, k, estKmercov2, estLength2, max_iterations)
```
