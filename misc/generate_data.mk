VARIABLES = genomescope_constants
DATA_FILES = $(patsubst %, data/%.rdata, $(VARIABLES))
# replace path by Rscript -e "noquote(.libPaths())" | tail -1 | cut -f 2 -d ' '
# INSTALLATION = /Library/Frameworks/R.framework/Versions/3.3/Resources/library/GenomeScope

.PHONY : prepare_data
prepare_data : $(DATA_FILES)

data/%.rdata : data-raw/%.R
	Rscript $<
