# alorenzetti 202012

# description #####
# this script will parse
# tpms generated by runKallisto
# pipeline

# defining function to compute se
rowSes = function(M){
  M = M %>% as.matrix()
  n = dim(M)[2]
  sds = M %>% rowSds()
  ses = sds / sqrt(n)
  
  return(ses)
}

# reading files ####
# exploratory analysis of kallisto TPM datasets
# to manually compute TPMs, follow the instructions on the following
# page: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# reading and parsing totalrna tpm (version to be used in the exploratory analysis)
totrnatpmraw=read_delim("data/tableTpmTotalRNA23_v3.tsv",delim="\t") %>% 
  select(-length,-eff_length) %>% 
  rename_with(.cols = -target_id,
              .fn = ~ str_replace(.x,
                                  "total-RNA-(.*)-(.*)_S.*$",
                                  "lysate_TP\\2_BR\\1"))

# version to be further used in all downstream analyses except exploratory analysis
totrnatpm=read_delim("data/tableTpmTotalRNA23_v3.tsv",delim="\t")
totrnatpm["mean_abundance_rna_total_TP1"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-1")) %>% rowMeans()
totrnatpm["mean_abundance_rna_total_TP2"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-2")) %>% rowMeans()
totrnatpm["mean_abundance_rna_total_TP3"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-3")) %>% rowMeans()
totrnatpm["mean_abundance_rna_total_TP4"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-4")) %>% rowMeans()

totrnatpm["se_abundance_rna_total_TP1"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-1")) %>% rowSes()
totrnatpm["se_abundance_rna_total_TP2"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-2")) %>% rowSes()
totrnatpm["se_abundance_rna_total_TP3"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-3")) %>% rowSes()
totrnatpm["se_abundance_rna_total_TP4"] = totrnatpm %>%
  dplyr::select(matches("total-RNA-[1-3]-4")) %>% rowSes()

totrnatpm = totrnatpm %>%
  dplyr::select(target_id, contains("abundance")) %>% 
  mutate(locus_tag = sub("\\|.*$", "", target_id)) %>% 
  dplyr::select(-target_id)

# adjusting colnames
totrnatpmlong = totrnatpm %>%
  dplyr::select(locus_tag, contains("mean")) %>% 
  pivot_longer(cols = contains("abundance"),
               names_to = c("measure", "libtype", "timepoint"),
               names_pattern = "^(.*)?_(.*)_(.*)$",
               values_to = "abundance")

# plotting abundance distribution for each timepoint
# ggplot(totrnatpmlong, aes(x=log10(abundance), color = timepoint)) +
#   geom_density() +
#   facet_grid(~ libtype)

# reading and parsing riboseq tpm (version to be used in the exploratory analysis)
ribornatpmraw = read_delim("data/tableTpmRiboSeqTrim15_v3.tsv",delim="\t") %>% 
  select(-length,-eff_length) %>% 
  rename_with(.cols = -target_id,
              .fn = ~ str_replace(.x,
                                  "ribosomal_RNA_(.*)-(.*)_S.*$",
                                  "ribo_TP\\2_BR\\1"))

# version to be further used in all downstream analyses except exploratory analysis
ribornatpm=read_delim("data/tableTpmRiboSeqTrim15_v3.tsv",delim="\t")
ribornatpm["mean_abundance_rna_ribofraction_TP1"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-1")) %>% rowMeans()
ribornatpm["mean_abundance_rna_ribofraction_TP2"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-2")) %>% rowMeans()
ribornatpm["mean_abundance_rna_ribofraction_TP3"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-3")) %>% rowMeans()
ribornatpm["mean_abundance_rna_ribofraction_TP4"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-4")) %>% rowMeans()

ribornatpm["se_abundance_rna_ribofraction_TP1"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-1")) %>% rowSes()
ribornatpm["se_abundance_rna_ribofraction_TP2"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-2")) %>% rowSes()
ribornatpm["se_abundance_rna_ribofraction_TP3"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-3")) %>% rowSes()
ribornatpm["se_abundance_rna_ribofraction_TP4"] = ribornatpm %>%
  dplyr::select(matches("ribosomal_RNA_[1-3]-4")) %>% rowSes()

ribornatpm = ribornatpm %>%
  dplyr::select(target_id, contains("abundance")) %>% 
  mutate(locus_tag = sub("\\|.*$", "", target_id)) %>% 
  dplyr::select(-target_id)

ribornatpmlong = ribornatpm %>%
  dplyr::select(locus_tag, contains("abundance")) %>% 
  pivot_longer(cols = contains("abundance"),
               names_to = c("measure", "timepoint"),
               names_pattern = "^(.*)?_.*_(.*)$",
               values_to = "abundance")

# plotting abundance distribution for each timepoint
# ggplot(ribornatpmlong, aes(x=log10(abundance), color = timepoint)) +
#   geom_density() +
#   facet_grid(~ measure)

# merging totrna, and riborna tpms
tpm = left_join(totrnatpm, ribornatpm, by = "locus_tag") %>% 
  relocate("locus_tag")

# following chunk is probably not needed anymore
# # regarding gvp1a cluster: I will change
# # the locus_tags to match those of gvp1b
# # gvp1a e gvp1b are copies
# # gvp2 = paste0("VNG_",
# #               c("6226a",
# #                 "6229G",
# #                 "6230G",
# #                 "6232G",
# #                 "6233G",
# #                 "6235G",
# #                 "6236G",
# #                 "6237G",
# #                 "6239G",
# #                 "6240G",
# #                 "6241G",
# #                 "6242G",
# #                 "6244G",
# #                 "6246G"))
# 
# gvp1a = paste0("VNG_",
#                c("7015",
#                  "7016",
#                  "7017",
#                  "7018",
#                  "7019",
#                  "7020",
#                  "7021",
#                  "7022",
#                  "7023",
#                  "7024",
#                  "7025",
#                  "7026",
#                  "7027",
#                  "7028"))
# 
# gvp1b = paste0("VNG_",
#                c("6019G",
#                  "6020G",
#                  "6021G",
#                  "6022G",
#                  "6023G",
#                  "6024G",
#                  "6025G",
#                  "6026G",
#                  "6027G",
#                  "6028G",
#                  "6029G",
#                  "6031G",
#                  "6032G",
#                  "6033G"))
# 
# gvp1bnames = paste0("(Gvp",
#                     c(LETTERS[13:4],
#                     c("A", "C", "N", "O")),
#                       ")")
# 
# tpm$locus_tag[grepl(paste0("^", gvp1a, "$") %>% paste0(collapse = "|"), tpm$locus_tag)] = gvp1b
# 
# # a few names in the protein dataset don't match
# # those of the rna-seq dataset; those are going
# # to be adjusted manually
# oldnames = c("VNG_7001",
#              "VNG_7002",
#              "VNG_7006",
#              "VNG_7007",
#              "VNG_7009",
#              "VNG_7029",
#              "VNG_7030",
#              "VNG_7037c",
#              "VNG_7038a",
#              "VNG_7039",
#              "VNG_7041",
#              "VNG_7042",
#              "VNG_7043",
#              "VNG_7048",
#              "VNG_7054",
#              "VNG_7056",
#              "VNG_7057a",
#              "VNG_7060",
#              "VNG_7061",
#              "VNG_7062",
#              "VNG_7064",
#              "VNG_7065",
#              "VNG_7066",
#              "VNG_7067",
#              "VNG_7069",
#              "VNG_7072",
#              "VNG_7075",
#              "VNG_7081",
#              "VNG_7082",
#              "VNG_7083",
#              "VNG_7085a",
#              "VNG_7086",
#              "VNG_7087",
#              "VNG_7090",
#              "VNG_7091",
#              "VNG_7092",
#              "VNG_7093")
# 
# newnames = c("VNG_6001H",
#              "VNG_6003H",
#              "VNG_6008H",
#              "VNG_6009H",
#              "VNG_6011H",
#              "VNG_6034G",
#              "VNG_6035G",
#              "VNG_6047H",
#              "VNG_6051H",
#              "VNG_6053G",
#              "VNG_6057C",
#              "VNG_6059C",
#              "VNG_6060C",
#              "VNG_6065G",
#              "VNG_6073G",
#              "VNG_6076H",
#              "VNG_6077H",
#              "VNG_6081G",
#              "VNG_6082H",
#              "VNG_6083H",
#              "VNG_6085H",
#              "VNG_6086G",
#              "VNG_6087C",
#              "VNG_6088C",
#              "VNG_6090C",
#              "VNG_6095C",
#              "VNG_6101H",
#              "VNG_6113H",
#              "VNG_6115H",
#              "VNG_6116H",
#              "VNG_6120H",
#              "VNG_6121H",
#              "VNG_6123G",
#              "VNG_6127H",
#              "VNG_6128H",
#              "VNG_6129C",
#              "VNG_6130G")
# 
# tpm$locus_tag[grepl(paste0("^", oldnames, "$") %>% paste0(collapse = "|"), tpm$locus_tag)] = newnames
