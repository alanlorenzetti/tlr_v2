# alorenzetti 202008

# description ####
# this script will get/load
# a few additional misc features
# LSm binding sites
# UTR info
# GC content

# loading libraries #####
source("scripts/loadingLibs.R")

# processing starts ####
# getting and parsing LSm binding info
# here I will load the table, remove the asRNAs
# and lookup into genes and 5UTRs
# if LSm binding sites are detected
# in genes or 5UTRs, sense or antisense,
# that will be taken into consideration
lsmGenes = read_delim(file = "data/resultsLSmGenes.csv", delim = ",") %>% 
  mutate(name = sub("^5UTR_", "", name)) %>% 
  select(-chr,-start,-end,-strand) %>% 
  group_by(name) %>% 
  summarise_all(.funs = list(sum)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(lsmSense = sum(one_sense, two_sense, three_sense, four_sense, five_sense),
         lsmAntiSense = sum(one_as, two_as, three_as, four_as, five_as)) %>% 
  ungroup() %>% 
  mutate(lsmSense = case_when(lsmSense > 0 ~ "yes",
                              TRUE ~ "no"),
         lsmAntiSense = case_when(lsmAntiSense > 0 ~ "yes",
                                  TRUE ~ "no")) %>% 
  dplyr::select(name, lsmSense, lsmAntiSense)

# getting and parsing info about UTRs
utr = read_delim(file = "data/5UTR.txt",
                 delim = "\t",
                 col_names = T) %>% 
  select(-length)

# getting and parsing info about UTR mfe
utrmfe = read_delim(file = "data/5UTRplus18ntdownstream_mfe.txt",
                    delim = "\t",
                    col_names = T) %>% 
  select(-seq,-dotbracket)

# getting GC content info for each gene
# based on sequences stored in pfeiSeq object
GCcontent = tibble(locus_tag = names(pfeiSeqs) %>%
                     sub("\\|.*$", "", .),
                   GC = pfeiSeqs %>%
                     letterFrequency(letters = "GC", as.prob = T) %>%
                     as.numeric(),
                   GCdev = GC - mean(GC))
