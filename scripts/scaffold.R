# alorenzetti 202008

# this is a scaffold script
# here you can find the adequate
# script order

# loading libs
source("scripts/loadingLibs.R")

# ~/gdrive/runKallisto/plotKmerTests.R
source("scripts/locusTagDictGbk2Tpa.R")
source("scripts/parseHalfLives.R")
source("scripts/codonUsage.R")
source("scripts/parseMicrobesOnlinePredOperons.R")
source("scripts/parseAntisenseRNAs.R")
source("scripts/getAndParseKEGG.R")
source("scripts/getAndParseUniprotInfo.R")
source("scripts/getAndParseArCOGs.R")
source("scripts/getAndParseMiscFeatures.R")
source("scripts/combineAndWrangleFunCat.R")

source("scripts/parseProteomicsSpectronaut.R")
source("scripts/parseProteomicsOneOmics.R")

source("scripts/parseRNAtpms.R")
source("scripts/DEanalysis.R")

source("scripts/models.R")
source("scripts/proteinRegRules.R")
source("scripts/timecourseAnalysis.R")

source("scripts/heatmaps.R")
source("scripts/clusterEnrichmentAnalysis.R")
source("scripts/unifyAbundance.R")
source("scripts/analysesAndFigures.R")