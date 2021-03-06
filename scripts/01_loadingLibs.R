# alorenzetti 202007

# description #####
# this script will load all the required packages
# for the downstream steps
# it will load essential files and objx as well
if(!require(pacman)){install.packages("pacman")}
library(pacman)

# this one isnt stored neither in cran nor in bioconductor
# https://github.com/drostlab/orthologr
# devtools::install_github("https://github.com/drostlab/orthologr", build_vignettes = T)
library(orthologr)

# the following is not available at CRAN
# library(devtools) ; install_github("allydunham/tblhelpr")
library(tblhelpr)

# required packages
packs = c("BiocManager",
          "tidyverse",
          "ComplexHeatmap",
          "circlize",
          "rtracklayer",
          "Biostrings",
          "GenomicRanges",
          "GenomicFeatures",
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
          "ggridges",
          "openxlsx",
          "mclust",
          "purrr",
          "svglite",
          "ggtext",
          "gtools",
          "DT",
          "clipr")

# loading and installing missing packages
p_load(char = packs)

# setting gglot2 theme
theme_set(theme_bw())

# loading tab10 palette
tab10 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`
tab20 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 20`
seqblues = ggthemes_data$tableau$`color-palettes`$`ordered-sequential`$Blue$value

# I have standardized a non-redundant transcriptome for
# Hsalinarum NRC-1. It is available at:
# https://alanlorenzetti.github.io/halo_nr_tx/
# downloading the list of non-redundant locus_tags
# with alternative names
nrtx = read_tsv("https://alanlorenzetti.github.io/halo_nr_tx/data/dictionary.tsv")
nrtxseqs = readDNAStringSet(filepath = "https://alanlorenzetti.github.io/halo_nr_tx/data/Hsalinarum_nrtx.fa",
                            format = "fasta")
