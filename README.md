# GenomeTelescope






### Thoughts for the development of the package

**models as classes**

It would be really handy if each of the models we will include would be a class, where each of them would their own functions:
   - predict 1n
	 - predict 2n
	 - predict
	 - plot

They still should stay as a subclass of the model, so all the higher level functions still work on it (confint, coef, etc)

**input of models**

I think the models take too many parameters and it's getting clunky.

We could also define a class for the k-mer histogram. The k-mer histogram would probably have just list with a `data.frame` (or two vectors) for frequency and coverage and some useful metadata, perhaps including starting values? Then a common check on each model could be if the starting value of the fitted parameter is or is not defined. Perhpas there could be two inner lists with fixed and fitted parameters (to make a clear distinction between them). So, the input of the model would be only the "kmerHist" class (or however we would call it).

Also, maybe we should avoid the word "fixed", because it's nothing like the "fixed" effects in mixed models.
