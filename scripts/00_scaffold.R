# alorenzetti 202008

# description ####
# this is a scaffold script
# here you can find the adequate
# script order

# creating a few directories necessary to store results and plots
if(!dir.exists("plots")){dir.create("plots")}
if(!dir.exists("results")){dir.create("results")}

# loading libs
source("scripts/01_loadingLibs.R")
source("scripts/02_parseHalfLives.R")
source("scripts/03_codonUsage.R")
source("scripts/04_parseMicrobesOnlinePredOperons.R")
source("scripts/05_parseAntisenseRNAs.R")
source("scripts/06_parseTFChipSeq.R")
source("scripts/07_parseTFChipChip.R")
source("scripts/08_getAndParseMiscFeatures.R")
source("scripts/09_combineAndWrangleFunCat.R")
source("scripts/10_parseProteomicsSpectronaut.R")
source("scripts/11_parseRNAtpms.R")
source("scripts/12_exploratoryAnalysis.R") # exploratory analysis plots for PhD thesis
source("scripts/13_unifyAbundance.R")
source("scripts/14_heatmapsAllAbsolute_v2_tese.R") # heat maps for PhD thesis
source("scripts/15_analysisAndFiguresRequested20201009.R") # correlation plots for PhD thesis