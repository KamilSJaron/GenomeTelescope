# GenomeTelescope

This is a package that should eventually allow us to do genome modeling in a more versatile fashion. It is fully based on the idea of GenomeScope, but it explores a somewhat broader range of scenarios where these models can be utilised.

### Background & motivation

Genomescope is awesome, but lately I run into problems where GenomeScope was unable to do what I wanted to. Why? Well, because not all sequenced genomes follow the assumption of GenomScope.

This repo will deal with very specific problem, if you just want a genome model, I would recommend to look at the [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0) or the original [GenomeScope](https://github.com/schatzlab/genomescope) instead.

This a split-fork of the older GenomScope version for one reason - it's a lot simpler code to tweak, but I will certainly use some of the tricks from GenomeScope 2.0 too (and possibly completely switch to a GenomeScope 2.0 fork). Also, this is a construction zone, I will be pushing commits with my progress, not everything will be necesarily working.

#### References

- Vurture, GW, Sedlazeck, FJ, Nattestad, M, Underwood, CJ, Fang, H, Gurtowski, J, Schatz, MC (2017) *Bioinformatics* doi: https://doi.org/10.1093/bioinformatics/btx153
- Ranallo-Benavidez TR, Jaron KS, Schatz MC. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. Nature Communications. 2020 Mar 18;11(1):1-0.

#### List of problems to be explored using multiple k-mer spectra

- uneven spacing of k-mer peaks
- the genome has a haploid portion (e.g. genome of the heterogametic sex)
   - X0 systems (easish)
   - XY/ZW systems (harder to converge)
   - system agnostic model to est differences in genome sizes (differences in genome sizes between libraries)


## Development notes and thoughts

**Testing datasets during development**

We need to figure out the best way to have a battery of testing k-mer histograms.

For now i created a makefile that can fetch us the data from pre-defined urls. To dl the files, run (expecting `wget` to be installed at you OS)

```
make -f misc/download_testing_data.mk
```

will download a bunch of kmer histograms in `data` subdirectory (which is in `.gitignore` so you don't have to be worried about accidently commiting them).

Of course, for smarter testing we should figure out how to include this data in the automated tests of the package (or alternativelly we can have a make-based integrative testing, but I would not be a big fan of that).

**consideration of multiple samples**

There is no good reason why NOT to include more samples if they are available. However, it is hard to figure out how we can actually deploy such information. E.g. imagine we have 4 males and 4 females and we are interested in knowing size of X and Y. Is there anything smarter than a bunch of pairwise comparisons?

**models as classes**

It would be really handy if each of the models we will include would be a class, where each of them would their own functions:
	 - predict (with argument specifying full, monosomic, or disomic portion)
	 - plot
   - print (so users have very clear access to the models excetly as they are fit)

They still should stay as a subclass of the model, so all the higher level functions still work on it (confint, coef, etc)

**input of models**

I think the models take too many parameters and it's getting clunky.

We could also define a class for the k-mer histogram. The k-mer histogram would probably have just list with a `data.frame` (or two vectors) for frequency and coverage and some useful metadata, perhaps including starting values? Then a common check on each model could be if the starting value of the fitted parameter is or is not defined. Perhpas there could be two inner lists with fixed and fitted parameters (to make a clear distinction between them). So, the input of the model would be only the "kmerHist" class (or however we would call it).

Also, maybe we should avoid the word "fixed", because it's nothing like the "fixed" effects in mixed models.
