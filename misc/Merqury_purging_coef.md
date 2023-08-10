K-mer histogram files

```bash
ilTorViri5.ccs.ilTorViri5.a_ctg.spectra-cn.cni
ilTorViri5.ccs.ilTorViri5.p_ctg.spectra-cn.cni
```

Let's try to fit a simple model

```R
library('GenomeTelescope')

ilTorViri5a <- read.table('data/Tvir/ilTorViri5.ccs.ilTorViri5.a_ctg.spectra-cn.cni', header = TRUE)
head(ilTorViri5a)

purged_hist <- ilTorViri5a[ilTorViri5a$Copies == "read-only", c(2, 3)]
single_copy_hist <- ilTorViri5a[ilTorViri5a$Copies == "1", c(2, 3)]

cov_range <- 15:100
coverage_barplot(purged_hist[cov_range, 2], purged_hist[cov_range, 1])
coverage_barplot(single_copy_hist[cov_range, 2], single_copy_hist[cov_range, 1])

x_purged <- purged_hist[cov_range, 1]
x_single <- single_copy_hist[cov_range, 1]
y_purged <- purged_hist[cov_range, 2]
y_single <- single_copy_hist[cov_range, 2]

estKmercov <- 30
estLength <- 400e6
max_iterations <- 20

# basic model of sorts - genomescope style
sinlge_copy_model <- nls(y_single ~ (a       * dnbinom(x_single, size = kmercov   / bias, mu = kmercov) +
                                     (1 - a) * dnbinom(x_single, size = kmercov*2 / bias, mu = kmercov * 2)) * length1,
                         start = list(a=0.2, kmercov=estKmercov, bias = 0.5, length1=estLength),
                         control = list(minFactor=1e-12, maxiter=20))

purged_model <- nls(y_purged ~ dnbinom(x_purged, size = kmercov   / bias, mu = kmercov) * length2,
                         start = list(kmercov=estKmercov, bias = 0.5, length2=estLength),
                         control = list(minFactor=1e-12, maxiter=20))

single_copy_monoploid = coef(sinlge_copy_model)['a'] * coef(sinlge_copy_model)['length1']
purged_monoploid = coef(purged_model)['length2']

single_copy_monoploid / purged_monoploid
```

TODO:
- fit just one model that would force the same k-mer coverage, and overdispersal in both
- consider also the 2n