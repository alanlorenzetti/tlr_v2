# alorenzetti 202007

# description #####
# this script will take deAlmeida et al. 2019
# annotation of antisense RNAs and find what
# genes of pfeiffer et al. 2019 intersect
# with them. this is part of functional
# categorization step

# loading libs ####
source("scripts/loadingLibs.R")

# reading antisense gff
asrnas = rtracklayer::import("data/Hsalinarum-846asRNAs-deAlmeida2019.gff3")
genes = rtracklayer::import("data/Hsalinarum-gene-annotation-pfeiffer2019.gff3")
genes = genes[genes$type == "gene",]

# inverting asrnas strand to find intersections
strand(asrnas) = invertStrand(strand(asrnas))

# checking intersections
geneswasrnas = subsetByOverlaps(genes, asrnas, ignore.strand=F) %>%
  as_tibble() %>%
  dplyr::select(locus_tag) %>%
  mutate(asRNA = "yes")
