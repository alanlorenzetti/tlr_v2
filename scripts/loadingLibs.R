# alorenzetti 202007

# description #####
# this script will load all the required packages
# for the downstream steps
if(!require(pacman)){install.packages("pacman")}
library(pacman)

# this one isnt stored neither in cran nor in bioconductor
# https://github.com/drostlab/orthologr
# devtools::install_github("https://github.com/drostlab/orthologr")
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
