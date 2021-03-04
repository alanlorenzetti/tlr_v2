# alorenzetti 202008

# this is a scaffold script
# here you can find the adequate
# script order

# loading libs
source("scripts/loadingLibs.R")

# ~/gdrive/runKallisto/plotKmerTests.R
#source("scripts/locusTagDictGbk2Tpa.R") # about to get deprecated; replaced by https://alanlorenzetti.github.io/halo_nr_tx/
source("scripts/parseHalfLives.R")
source("scripts/codonUsage.R")
source("scripts/parseMicrobesOnlinePredOperons.R")
source("scripts/parseAntisenseRNAs.R")
#source("scripts/getAndParseKEGG.R") # about to get deprecated
#source("scripts/getAndParseUniprotInfo.R") # about to get deprecated
#source("scripts/getAndParseArCOGs.R") # about to get deprecated
source("scripts/parseTFChipSeq.R")
source("scripts/parseTFChipChip.R")
source("scripts/getAndParseMiscFeatures.R")
source("scripts/combineAndWrangleFunCat.R")

source("scripts/parseProteomicsSpectronaut.R")
# source("scripts/parseProteomicsOneOmics.R") # about to get deprecated

source("scripts/parseRNAtpms.R")
# source("scripts/DEanalysis.R") # about to get deprecated

source("scripts/exploratoryAnalysis.R") # exploratory analysis plots for dissertation

# source("scripts/models.R") # about to get deprecated
# source("scripts/proteinRegRules.R") # about to get deprecated
# source("scripts/timecourseAnalysis.R") # about to get deprecated

#source("scripts/heatmaps.R")
# source("scripts/clusterEnrichmentAnalysis.R") # about to get deprecated
source("scripts/unifyAbundance.R")
# source("scripts/heatmapsAllRelativeChanges.R") # about to get deprecated
# source("scripts/heatmapsAllAbsolute.R") # about to get deprecated
source("scripts/heatmapsAllAbsolute_v2_tese.R")
# source("scripts/computeAnglesAndMag.R") # about to get deprecated
source("scripts/analysesAndFigures.R")
source("scripts/analysisAndFiguresRequested20201009.R") # correlation plots
source("scripts/findPutativeRegulatedTx.R")
# source("scripts/norm_issues_tmp.R") # accessory and not essential
# source("scripts/bepeReportFigures.R")