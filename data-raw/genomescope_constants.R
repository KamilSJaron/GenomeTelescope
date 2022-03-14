## Number of rounds before giving up
NUM_ROUNDS=4

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Max rounds on NLS
MAX_ITERATIONS=20

## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
SCORE_CLOSE = 0.20

## Overrule heterozygosity if there is a large difference in het rate
SCORE_HET_FOLD_DIFFERENCE = 10

## Print out VERBOSEging messages (0/1)
VERBOSE = 0

## Colors for plots
COLOR_BGCOLOR  = "light grey"
COLOR_HIST     = "#56B4E9"
COLOR_4PEAK    = "black"
COLOR_2PEAK    = "#F0E442"
COLOR_ERRORS   = "#D55E00"
COLOR_KMERPEAK = "black"
COLOR_RESIDUAL = "purple"
COLOR_COVTHRES = "red"

save(list = ls(), file = "data/genomescope_constants.RData")
