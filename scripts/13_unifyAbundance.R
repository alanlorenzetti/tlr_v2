# alorenzetti 202008

# description #####
# this script is intended to build
# linear models and figures
# using abundance data

# unifying protein counts and rna counts ####
# and performing quantile normalization
abund = left_join(spectroWide, tpm, by = "locus_tag")

abundmean = abund %>%
  dplyr::select(starts_with("mean")) %>%
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble(rownames = NULL)
colnames(abundmean) = abund %>%
  dplyr::select(starts_with("mean")) %>% 
  colnames()

abundse = abund %>%
  dplyr::select(starts_with("se")) %>%
  as.matrix() %>% 
  normalize.quantiles() %>% 
  as_tibble(rownames = NULL)
colnames(abundse) = abund %>%
  dplyr::select(starts_with("se")) %>% 
  colnames()

abundNorm = bind_cols(abund[,"locus_tag"],
                      abundmean,
                      abundse)

abundNormDer = bind_cols(abund[,"locus_tag"],
                         abundmean,
                         abundse) %>% 
  mutate(mean_abundance_rna_occupancy_TP1 = mean_abundance_rna_ribofraction_TP1 / mean_abundance_rna_total_TP1,
         mean_abundance_rna_occupancy_TP2 = mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_total_TP2,
         mean_abundance_rna_occupancy_TP3 = mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_total_TP3,
         mean_abundance_rna_occupancy_TP4 = mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_total_TP4) %>%
  
  mutate(mean_abundance_rna_psiTE_TP1 = mean_abundance_protein_lysate_TP1 / mean_abundance_rna_total_TP1,
         mean_abundance_rna_psiTE_TP2 = mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP2,
         mean_abundance_rna_psiTE_TP3 = mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP3,
         mean_abundance_rna_psiTE_TP4 = mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP4) %>%
  
  mutate(mean_abundance_rna_tlr_TP1 = mean_abundance_rna_psiTE_TP1 / mean_abundance_rna_occupancy_TP1,
         mean_abundance_rna_tlr_TP2 = mean_abundance_rna_psiTE_TP2 / mean_abundance_rna_occupancy_TP2,
         mean_abundance_rna_tlr_TP3 = mean_abundance_rna_psiTE_TP3 / mean_abundance_rna_occupancy_TP3,
         mean_abundance_rna_tlr_TP4 = mean_abundance_rna_psiTE_TP4 / mean_abundance_rna_occupancy_TP4)

# pivoting a version of abund non normalized dataset
abundLong = abund %>% 
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# pivoting a version of abund normalized dataset
abundNormLong = abundNorm %>% 
  mutate(mean_abundance_protein_lysate_TP0 = 0,
         mean_abundance_protein_ribo_TP0 = 0,
         mean_abundance_rna_total_TP0 = 0,
         mean_abundance_rna_ribofraction_TP0 = 0) %>%
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# pivoting a version of abund normalized dataset with derived vars
abundNormDerLong = abundNormDer %>% 
  mutate(mean_abundance_rna_occupancy_TP0 = 1,
         mean_abundance_rna_psiTE_TP0 = 1,
         mean_abundance_rna_tlr_TP0 = 1) %>%
  pivot_longer(cols = contains("abundance"),
               names_pattern = "^(.*)_.*_(.*_.*)_(.*)$",
               names_to = c("measure", "libtype", "timepoint"),
               values_to = "abundance") %>% 
  pivot_wider(names_from = c("measure"),
              values_from = "abundance")

# # plots #####
# # densities
# abundLong %>%
#   ggplot(aes(x=log10(mean+1), color = libtype, linetype = timepoint)) +
#   geom_density()
# 
# abundNormLong %>%
#   ggplot(aes(x=log10(mean+1), color = libtype, linetype = timepoint)) +
#   geom_density()

# # trajectories
# abundLong %>%
# #  filter(!(libtype == "rna_psiTE" | libtype == "rna_occupancy")) %>% 
#   filter(locus_tag %in% "VNG_1496G") %>%
#   ggplot(aes(x=timepoint, y=mean, colour = libtype, group = libtype)) +
#   geom_line(size = 1, alpha = 0.75) +
#   geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
# #  facet_wrap(~locus_tag, "free_y") +
#   xlab("Time Point") + ylab("Abundance") +
#   scale_color_manual(values = c("rna_total" = "#E15759",
#                                 "rna_ribofraction" = "#F28E2B",
#                                 "protein_lysate" = "#4E79A7",
#                                 "protein_ribo" = "#59A14F"),
#                      name = "Variable")
# 
# abundNormLong %>%
# #  filter(!(libtype == "rna_psiTE" | libtype == "rna_occupancy")) %>% 
#   filter(locus_tag %in% "VNG_1496G") %>%
#   ggplot(aes(x=timepoint, y=mean, colour = libtype, group = libtype)) +
#   geom_line(size = 1, alpha = 0.75) +
#   geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
# #  facet_wrap(~locus_tag, "free_y") +
#   xlab("Time Point") + ylab("Abundance") +
#   scale_color_manual(values = c("rna_total" = "#E15759",
#                                 "rna_ribofraction" = "#F28E2B",
#                                 "protein_lysate" = "#4E79A7",
#                                 "protein_ribo" = "#59A14F"),
#                      name = "Variable")

# heatmap
# colfunct = circlize::colorRamp2(breaks = c(0,3,6), colors = viridis(3))
# M = abundNorm %>%
#   dplyr::select(-locus_tag,-starts_with("se")) %>%
#   dplyr::select(ends_with("TP1")) %>%
#   drop_na() %>%
#   as.matrix()
# M = M + 1
# M = log10(M)
# 
# Heatmap(M, col = colfunct)
