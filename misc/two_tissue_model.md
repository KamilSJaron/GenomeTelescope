## Set up

Original GenomeScope fits a model that is unsuitable for estimating 1n and 2n coverage indipendently. Which makes sense for most of the genomes. But globular springtail males do not have this luxury. Their body is formed of two tissues with different karyotypes and that is causing the spacing of the 1n and 2n deviate a lot from 1:2 ratio. I did a very thourough exploration of this dataset in particular, you can check [our perprint](https://www.biorxiv.org/content/10.1101/2021.11.12.468426v1). I estimated the position of 1n and 2n peaks in two male samples using reads mapped back to the assembly, but I think it should be possible to get this model fit directly to k-mer spectra, because that's the observation we initially made - the k-mmer spectra was off.

For this functionality, I will use springtail testing data that can be obtained by

```
make -f misc/download_testing_data.mk
```

## Analysis

```{R}
library('GenomeTelescope')
# read kmer hist

input_file1 <- 'data/springtails/Afus1_k21_truncated.hist' # high coverage male
input_file2 <- 'data/springtails/BH3-2_k21_truncated.hist' # low coverage male
input_file3 <- 'data/springtails/Ocin2_k21_truncated.hist' # high coverage male of species without strange patters
input_file4 <- 'data/springtails/WW5-3_k21_truncated.hist' # low coverage female

# <- nlsLM_2peak_two_tissue_model()
```
