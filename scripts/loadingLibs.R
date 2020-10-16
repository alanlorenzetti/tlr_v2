# alorenzetti 202007

# description #####
# this script will load all the required packages
# for the downstream steps
if(!require(pacman)){install.packages("pacman")}
library(pacman)

# this one isnt stored neither in cran nor in bioconductor
# https://github.com/drostlab/orthologr
# devtools::install_github("https://github.com/drostlab/orthologr", build_vignettes = T)
library(orthologr)

# required packages
packs = c("BiocManager",
          "tidyverse",
          "ComplexHeatmap",
          "circlize",
          "rtracklayer",
          "Biostrings",
          "GenomicRanges",
          "NbClust",
          "BSgenome",
          "ggbeeswarm",
          "coRdon",
          "UniProt.ws",
          "gage",
          "readxl",
          "DESeq2",
          "tximport",
          "htmlwidgets",
          "plotly",
          "preprocessCore",
          "eulerr",
          "ggpubr",
          "DescTools",
          "see",
          "grid",
          "gridtext",
          "ggthemes",
          "TMixClust",
          "viridis",
          "matrixStats",
          "scales",
          "ggridges")

# loading and installing missing packages
p_load(char = packs)

# setting gglot2 theme
theme_set(theme_bw())

# loading tab10 palette
tab10 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`
seqblues = ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Blue$value
