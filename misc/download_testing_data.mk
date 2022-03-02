# makefile for the testing data

# this is how we specify variables with multiple elements
SPRINGTAIL_SAMPLES := Afus1 BH3-2 Ocin2 WW5-3
# this is how we can expand the variable according to the specified pattern
SPRINGTAIL_FILES := $(patsubst %, data/springtails/%_k21_truncated.hist, $(SPRINGTAIL_SAMPLES))

## download_all: a master call for downloading everything
.PHONY: download_all
download_all: $(SPRINGTAIL_FILES)

# this is just to make sure the springtail directory exist
data/springtails:
	mkdir -p $@

# this is the part that downloads the actual histograms from a different (more spoiled) github repo
# the "data/springtails" dependency is the directory
data/springtails/%_k21_truncated.hist: data/springtails
	wget "https://raw.githubusercontent.com/KamilSJaron/genomescope/master/analysis/real_data/springtails/"$*"_k21_truncated.hist" -O $@
