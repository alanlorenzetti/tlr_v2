# alorenzetti 202008

# descripition ####
# in this script we will have
# figures and additional analysis

# loading libs #####
source("scripts/loadingLibs.R")

# adding functional categorization to
# abundance dataframes
abundLongFuncat = left_join(abundLong,
                            dictFunCat,
                            by = c("locus_tag" = "pfeiLocusTag"))

abundNormFuncat = left_join(abundNorm,
                            dictFunCat,
                            by = c("locus_tag" = "pfeiLocusTag"))

# Lsm trajectory
abundLongFuncat %>% 
  filter(locus_tag %in% "VNG_1496G") %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# ISH2 trajectory
abundLongFuncat %>% 
  filter(locus_tag %in% "VNG_0210H") %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# cluster proteinDown mRNA Up in the context of LSm protein
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=regRule, fill = lsmSense)) + 
  geom_bar(position = "fill") +
  coord_flip()

# cluster proteinDown mRNA Up in the context of total mRNA abundance
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = regRule, x = log10(totrnaBaseMean))) +
  geom_density_ridges()

# lsm in context of HL
abundLongFuncat %>% 
  filter(libtype == "rna_total") %>% 
  ggplot(aes(x = log10(mean), y = HL, color = lsmSense)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  xlab("log10(RNA Abundance)")

# regRule in context of HL
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = regRule, x = HL)) +
  geom_density_ridges()

# regRule in context of cai
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = regRule, x = cai)) +
  geom_density_ridges()
  
# cluster proteinDown mRNA Up trajectory
pDmUTP4 = joinedTibble[["TP4_vs_TP4"]] %>% 
  filter(regRule == "ProteinDown_mRNAUp") %>% 
  select(locus_tag) %>% 
  unlist(use.names = F)

abundLongFuncat %>% 
  filter(locus_tag %in% pDmUTP4) %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# cluster mobilome in the context of LSm protein
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=arCOG, fill = lsmSense)) + 
  geom_bar(position = "fill") +
  coord_flip()

# cluster mobilome in the context of total mRNA abundance
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(y = arCOG, x = log10(totrnaBaseMean))) +
  geom_density_ridges()

# cluster mobilome in the context of asRNA
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=arCOG, fill = asRNA)) + 
  geom_bar(position = "fill") +
  coord_flip()

# cluster mobilome in the context of regRules
joinedTibble[["TP4_vs_TP4"]] %>% 
  ggplot(aes(x=arCOG, fill = regRule)) + 
  geom_bar(position = "fill") +
  coord_flip()

# lsm in context of totalmRNA abundance
abundLongFuncat %>% 
  filter(libtype == "rna_total") %>% 
  ggplot(aes(x = log10(mean), color = lsmSense)) +
  geom_density() +
  facet_wrap(~ timepoint)

# cluster proteinDown mRNA Up trajectory
mobilomeTP4 = joinedTibble[["TP4_vs_TP4"]] %>% 
  filter(str_detect(arCOG, "^Mobilome")) %>% 
  select(locus_tag) %>% 
  unlist(use.names = F)

abundLongFuncat %>% 
  filter(locus_tag %in% mobilomeTP4) %>% 
  filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# exploring abundance variables ####
# protein abundance in function of mRNA
# color == asRNA
abundLongFuncat %>% 
select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log10(protein_lysate),
             x = log10(rna_total),
             color = asRNA)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log10(Protein Abundance)") +
  xlab("log10(mRNA Abundance)") 

# protein abundance in function of mRNA
# color == LSm
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log10(protein_lysate),
             x = log10(rna_total),
             color = lsmSense)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log10(Protein Abundance)") +
  xlab("log10(mRNA Abundance)") 

# proteins in function rpfs
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log10(protein_lysate),
             x = log10(rna_ribofraction))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log10(Protein Abundance)") +
  xlab("log10(RPF Abundance)")

# what is ribosome occupancy? ####
# are they active translating ribosomes
# or stalled ribosomes?
# occupancy vs total
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_occupancy),
             x = log10(rna_total),
             color = cai)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Ribosome Occupancy)") +
  xlab("log10(mRNA Abundance)")  +
  scale_color_distiller(name = "cai",
                        type = "div",
                        palette = "RdYlGn")

# occupancy vs protein abundance
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_occupancy),
             x = log10(protein_lysate))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Ribosome Occupancy)") +
  xlab("log10(Protein Abundance)")

# occupancy vs cai
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_occupancy),
             x = cai)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Ribosome Occupancy)") +
  xlab("CAI")

# what means translational efficiency? ####
# translational efficiency in function of total mRNA
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log10(rna_total))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log10(mRNA Abundance)")

# translational efficiency in function of protein abundance
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log10(protein_lysate))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log10(Protein Abundance)")

# translational efficiency in function of rpf abundance
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log10(rna_ribofraction))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log10(RPF Abundance)")

# translational efficiency in function of occupancy
abundLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_psiTE),
             x = log2(rna_occupancy))) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("log2(Protein / mRNA)") +
  xlab("log2(Ribosome Occupancy)")
  
# checking if abundance and LSm are related
abundLongFuncat %>% 
  ggplot(aes(y=log10(mean), x=timepoint, fill = libtype)) +
  geom_violin() +
  facet_grid(~ lsmSense)

# checking if abundance and asRNA are related
abundLongFuncat %>% 
  ggplot(aes(y=log10(mean), x=timepoint, fill = libtype)) +
  geom_boxplot() +
  facet_grid(~ asRNA)

# combining both with facets
abundLongFuncat %>% 
  ggplot(aes(y=log10(mean), x=timepoint, fill = libtype)) +
  geom_boxplot() +
  facet_grid(lsmSense ~ asRNA)

# checking if abundance and cai are related
abundLongFuncat %>% 
  ggplot(aes(x = cai, y = log10(mean), color = timepoint)) +
  geom_point() +
  facet_grid(~ libtype)

# checking trajectories of ribosome proteins
abundLongFuncat %>% 
  filter(locus_tag %in% riboProts) %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# checking trajectories of gas vesicles
gasVesProts = dictFunCat$pfeiLocusTag[dictFunCat$pfeiProduct %>% str_detect("vesicle")]
abundLongFuncat %>% 
  filter(locus_tag %in% gasVesProts) %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# checking trajectories of mobilome proteins
abundLongFuncat %>% 
  filter(arCOG == "Mobilome: prophages, transposons") %>% 
  ggplot(aes(x = timepoint,
             y = log10(mean),
             color = libtype,
             group = libtype)) +
  geom_line() +
  geom_linerange(aes(ymin = log10(mean-se), ymax = log10(mean+se))) +
  facet_wrap(~ locus_tag)

# checking relationship between occupancy and protein abundance
# seems like the higher the occupancy, the lower is the protein content
# is RPFs reflecting stale ribosomes?
# that makes sense, considering no drug was added to stop the ribosomes
abundNormFuncat %>% 
  select(-starts_with("se")) %>% 
  ggplot(aes(x=log10(mean_abundance_protein_lysate_TP1), y=log2(mean_abundance_rna_occupancy_TP1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ lsmSense + asRNA)

# higher occupancy also reflects in low codon adaptation index
# suggesting occupancy is actually a measure of ribosomes stalled
abundNormFuncat %>% 
  select(-starts_with("se")) %>% 
  ggplot(aes(x=cai, y=log2(mean_abundance_rna_occupancy_TP4))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ lsmSense + asRNA)

# occupancy is also supposed to play
# a role in half lives
abundNormFuncat %>% 
  select(-starts_with("se")) %>% 
  ggplot(aes(x=HL, y=log2(mean_abundance_rna_occupancy_TP1))) +
  geom_point() +
  geom_smooth(method = "lm")

# what about GC content?
abundNormFuncat %>% 
  select(-starts_with("se")) %>% 
  ggplot(aes(x=GC, y=log2(mean_abundance_rna_occupancy_TP1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ lsmSense)

# checking what are the most abundant proteins in the lysate fraction
# that is going to show us a reference
# set that can be used to compute the codon adaptation index
mostAbundProtsLysate = abundNorm %>% 
  select(locus_tag, starts_with("mean_abundance_protein_lysate")) %>% 
  mutate(meanAcrossTP = rowMeans(.[,-1], na.rm = T)) %>% 
  select(locus_tag, meanAcrossTP) %>% 
  mutate(meanAcrossTP = case_when(is.nan(meanAcrossTP) ~ NA_real_,
                                  TRUE ~ as.numeric(meanAcrossTP)))

fil95 = mostAbundProtsLysate$meanAcrossTP %>% quantile(probs = .95, na.rm = T)
mostAbundProtsLysate = mostAbundProtsLysate$locus_tag[mostAbundProtsLysate$meanAcrossTP > fil95]
mostAbundProtsLysate = mostAbundProtsLysate[!is.na(mostAbundProtsLysate)]
dictFunCat[dictFunCat$pfeiLocusTag %in% mostAbundProtsLysate,c("pfeiLocusTag", "pfeiProduct", "arCOGproduct")] %>% view

# what are the most abundant proteins in the ribosome fraction proteome?
mostAbundProtsRibo = abundNorm %>% 
  select(locus_tag, starts_with("mean_abundance_protein_ribo")) %>% 
  mutate(meanAcrossTP = rowMeans(.[,-1], na.rm = T)) %>% 
  select(locus_tag, meanAcrossTP) %>% 
  mutate(meanAcrossTP = case_when(is.nan(meanAcrossTP) ~ NA_real_,
                                  TRUE ~ as.numeric(meanAcrossTP)))

fil95 = mostAbundProtsRibo$meanAcrossTP %>% quantile(probs = .95, na.rm = T)
mostAbundProtsRibo = mostAbundProtsRibo$locus_tag[mostAbundProtsRibo$meanAcrossTP > fil95]
mostAbundProtsRibo = mostAbundProtsRibo[!is.na(mostAbundProtsRibo)]
dictFunCat[dictFunCat$pfeiLocusTag %in% mostAbundProtsRibo,c("pfeiLocusTag", "pfeiProduct", "arCOGproduct")] %>% view

# checking what are the shared proteins in
# the ribosome and lysate fractions
plot(euler(compact(list(lysate = mostAbundProtsLysate,
           ribo = mostAbundProtsRibo))),
     quantities = T)

# 18 proteins in ribo not contained in lysate
dictFunCat[dictFunCat$pfeiLocusTag %in% (mostAbundProtsRibo[!mostAbundProtsRibo %in% mostAbundProtsLysate]),]
# another question: from the set of most abundant in lysate
# how many are ribosome proteins?

# 22 proteins in lysate not contained in ribo
dictFunCat[dictFunCat$pfeiLocusTag %in% (mostAbundProtsLysate[!mostAbundProtsLysate %in% mostAbundProtsRibo]),]
# another question: from the set of most abundant in ribo
# how many are ribosome proteins?

