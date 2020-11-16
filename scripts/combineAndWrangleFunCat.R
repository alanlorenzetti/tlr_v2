# alorenzetti 202008

# description ####
# this script will take parsed dataframes
# generated previously and unify
# in a comprenhensive dataframe

# loading libraries
source("scripts/loadingLibs.R")

# adding functional characterization
# combining data from uniprot with arCOG
uniArCOG = left_join(uniprotHalDf, arcogHal, by = c("arCOG_ID" = "arCOG_ID")) %>% 
  group_by(locus_tag, UNIPROTKB, GO) %>% 
  summarise(arCOG_ID = paste0(arCOG_ID,collapse = ";"),
            arCOGcode = paste0(arCOGcode,collapse = ";"),
            arCOG = paste0(arCOG,collapse = ";"),
            arCOGproduct = paste0(arCOGproduct,collapse = ";")) %>% 
  ungroup() %>% 
  mutate(arCOGcode = sub("^NA;|;NA;|;NA$", "", arCOGcode),
         arCOG = sub("^NA;|;NA;|;NA$", "", arCOG))
  
# only three genes have more than one distinct arCOG:
# VNG0958G: T;N
# VNG1131G: D;O
# VNG1951G: V;E
# I am arbitrarily accepting only the first category for each gene
uniArCOG = uniArCOG %>%
  mutate(arCOGcode = sub(";.*$","",arCOGcode),
         arCOG = sub(";.*$","",arCOG),
         arCOGproduct = sub(";.*$","",arCOGproduct))

# combining with kegg data
uniArCOGKEGG = left_join(uniArCOG, keggSetHalDf, by = c("locus_tag" = "locus_tag"))

# operons
uniArCOGKEGGOper = left_join(uniArCOGKEGG, operons, by = "locus_tag")

# adding products and alternative locus_tag
dictFunCat = left_join(dictProd, uniArCOGKEGGOper, by = c("query_id" = "locus_tag"))

# changing colnames
colnames(dictFunCat)[colnames(dictFunCat) == "subject_id"] = "pfeiLocusTag"
colnames(dictFunCat)[colnames(dictFunCat) == "product.y"] = "pfeiProduct"
colnames(dictFunCat)[colnames(dictFunCat) == "query_id"] = "locus_tag"

######################
# at this moment, dictFunCat has non repetitive locus_tags
# and repetitive pfeiLocus_tags. from now on
# I will create a functional categorization
# containing only non redundant pfei locus tags
# if one want to check correspondence of pfei locus tags
# to locus_tags, check dictProd object
# I will also drop a few cols that I am not using right now
dictFunCat = dictFunCat %>% 
  dplyr::select(-UNIPROTKB,-GO,-KEGGpathway,-operon) %>% 
  arrange(pfeiLocusTag,arCOG_ID)

# I am keeping only the first entry of repetitive pfeiLocusTags
dictFunCat = dictFunCat[!dictFunCat$pfeiLocusTag %>% duplicated(),]

# adding lsm data
dictFunCat = left_join(dictFunCat, lsmGenes, by = c("pfeiLocusTag" = "name"))

# adding UTR info
dictFunCat = left_join(dictFunCat, utr, by = c("pfeiLocusTag" = "locus_tag"))

# adding asRNA info
dictFunCat = left_join(dictFunCat, geneswasrnas, by = c("pfeiLocusTag" = "locus_tag"))

# adding utr mfe info
dictFunCat = left_join(dictFunCat, utrmfe, by = c("pfeiLocusTag" = "locus_tag")) %>% 
  rename(utrGC = GC,
         utrmfe = mfe)

# adding GC content info
dictFunCat = left_join(dictFunCat, GCcontent, by = c("pfeiLocusTag" = "locus_tag"))

# adding halflives from Hundt et al. 2007
dictFunCat = left_join(dictFunCat, halfLives, by = c("pfeiLocusTag" = "locus_tag"))

# adding codon index usage
dictFunCat = left_join(dictFunCat, cai, by = c("pfeiLocusTag" = "locus_tag"))

# adding transcription factor information obtained
# from chip-seq 
dictFunCat = left_join(dictFunCat, chipSeqTFs, by = c("pfeiLocusTag" = "locus_tag"))

# and chip-chip experiments
dictFunCat = left_join(dictFunCat, chipChipTFs, by = c("pfeiLocusTag" = "locus_tag"))

# adjusting missing values for a few vars
dictFunCat = dictFunCat %>% 
  mutate(lsmSense = case_when(is.na(lsmSense) ~ "no",
                              TRUE ~ as.character(lsmSense)),
         lsmAntiSense = case_when(is.na(lsmAntiSense) ~ "no",
                                  TRUE ~ as.character(lsmAntiSense)),
         arCOG = case_when(is.na(arCOG) ~ "Function unknown",
                           arCOG == "NA" ~ "Function unknown",
                           TRUE ~ as.character(arCOG)),
         utrSize = case_when(is.na(utrSize) ~ "Unknown",
                             TRUE ~ as.character(utrSize)),
         asRNA = case_when(is.na(asRNA) ~ "no",
                           TRUE ~ as.character(asRNA))
)

# adjusting missing values for TF-related vars
dictFunCat = dictFunCat %>% 
  mutate(across(starts_with("ChIP"),
                .fns = ~ case_when(is.na(.x) ~ "no",
                                   TRUE ~ as.character(.x))))

# I wondered why ISH2 was categorized as a defense mechanism in arCOG
# apparently, ISH2's arCOG09385 is described as NikR gene
# Transcriptional regulator, CopG/Arc/MetJ family (DNA-binding and a metal-binding domains)
# I feel like rather changing that class to mobilome
# also, all transposases according to pfeiProduct will receive mobilome class
dictFunCat = dictFunCat %>%
  mutate(arCOG = case_when(grepl(x = dictFunCat$pfeiProduct, pattern = "(ISH2)") ~ "Mobilome: prophages, transposons",
                           TRUE ~ as.character(arCOG)),
         arCOGcode = case_when(grepl(x = dictFunCat$pfeiProduct, pattern = "(ISH2)") ~ "X",
                               TRUE ~ as.character(arCOGcode))) %>%
  mutate(arCOG = case_when(grepl(x = dictFunCat$pfeiProduct, pattern = "transposase") ~ "Mobilome: prophages, transposons",
                           TRUE ~ as.character(arCOG)),
         arCOGcode = case_when(grepl(x = dictFunCat$pfeiProduct, pattern = "transposase") ~ "X",
                               TRUE ~ as.character(arCOGcode)))


########## VISUALIZATION ###############
# funcat visualization
# lsm per arcog
# the results show lsm is proportionally
# more represented in mobilome and defense mechanism
# category; we have to be cautious, because there
# are repetitive elements within this categories

# # lsm per arcog
# ggplot(data = dictFunCat) +
#   geom_bar(aes(x=arCOG, fill=lsmSense), position="fill") +
#   coord_flip() +
#   ggtitle("LSm bound RNAs and arCOGs")
# 
# # antisense per arcog
# ggplot(data = dictFunCat) +
#   geom_bar(aes(x=arCOG, fill=asRNA), position="fill") +
#   coord_flip() +
#   ggtitle("Antisense RNAs and arCOGs")
# 
# # utrSize per arcog
# ggplot(data = dictFunCat) +
#   geom_bar(aes(x=arCOG, fill=utrSize), position="fill") +
#   coord_flip() +
#   ggtitle("5' UTR length and arCOGs")
# 
# # # utr mfe vs lsm
# # ggplot(data = dictFunCat, aes(y=mfe, x=lsmSense)) +
# #   geom_violin() +
# #   geom_beeswarm()
# 
# ggplot(data = dictFunCat, aes(x=mfe, color=lsmSense)) +
#   geom_density() +
#   ggtitle("5' UTR MFE distributions")
# 
# # # utr mfe vs GC
# # ggplot(data = dictFunCat, aes(y=GC, x=lsmSense)) +
# #   geom_boxplot() +
# #   geom_beeswarm()
# 
# ggplot(data = dictFunCat, aes(x=GC, color=lsmSense)) +
#   geom_density() +
#   ggtitle("5' UTR GC content distributions")
# 
# # checking half lives per arcog
# ggplot(data = dictFunCat, aes(y=HL, x=arCOG)) +
#   geom_boxplot() +
# #  geom_beeswarm() + 
#   coord_flip() +
#   ggtitle("Half Life distribution per category")
# 
# # comparing if distributions of GC are
# # different for lsm bound and unbound UTRs
# # we are going to compute empirical cumulative
# # distribution function for both of them and
# # check if they are similar
# # lsm bound
# lsmbound = dictFunCat[dictFunCat$lsmSense == "no","GC"] %>% drop_na() %>% unlist() %>% unname()
# lsmboundFun = ecdf(lsmbound)
# 
# lsmunbound = dictFunCat[dictFunCat$lsmSense == "yes","GC"] %>% drop_na() %>% unlist() %>% unname()
# lsmunboundFun = ecdf(lsmunbound)
# 
# # plotting
# plot(ecdf(lsmbound), xlim=range(c(lsmbound, lsmunbound)))
# plot(ecdf(lsmunbound), add=TRUE, lty="dashed")
# 
# # performing kolmogorov-smirnoff to compare
# # both functions
# # and we reject the null hypothesis
# ks.test(x = lsmbound, y = lsmunboundFun)
