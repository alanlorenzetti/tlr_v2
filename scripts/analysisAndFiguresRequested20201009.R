# 20201009
# request was made on 20201009

# description #####
# Nitin has requested a couple of charts to show the correlation between 
# transcripts, RPF, and protein levels. Ideally, we should be able to
# identify cases of transcriptional regulation, translational efficiency
# and translational regulation

# loading libs #####
source("scripts/loadingLibs.R")

# plotting ####
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

# protein abundance in function of mRNA ####
# no colorspace  
protvsmrna = abundNormLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks()

# saving
ggsave(filename = "plots/11_prot_vs_mrna_tpwise.png",
       plot = protvsmrna,
       units = "in",
       width = 6,
       height = 5)

# slide prot from TPi+1 mRNA from TPi
# TP2 vs TP1
slide1 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP2,
             x = mean_abundance_rna_total_TP1)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. mRNA TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_total_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. mRNA TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_total_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP4 vs. mRNA TP3")

# arranging plots
protvsmrnaslide = ggarrange(slide1, slide2, slide3,
                            ncol = 3, nrow = 1)

# saving
ggsave(filename = "plots/12_prot_vs_mrna_slides.png",
       plot = protvsmrnaslide,
       units = "in",
       width = 10,
       height = 3)

# Protein in function of RPF ####
# no colorspace  
proteinvsrpf = abundNormLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_ribofraction)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("RPF Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks()

# saving
ggsave(filename = "plots/21_prot_vs_rpf_tpwise.png",
       plot = proteinvsrpf,
       units = "in",
       width = 6,
       height = 5)

# slide prot from TPi+1 rpf from TPi
# TP2 vs TP1
slide1 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP2,
             x = mean_abundance_rna_ribofraction_TP1)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("RPF Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. RPF TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_ribofraction_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("RPF Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. RPF TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_ribofraction_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("RPF Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks() +
  ggtitle("Protein TP4 vs. RPF TP3")

# arranging plots
protvsrpfslide = ggarrange(slide1, slide2, slide3,
                           ncol = 3, nrow = 1)

# saving
ggsave(filename = "plots/22_prot_vs_rpf_slides.png",
       plot = protvsrpfslide,
       units = "in",
       width = 10,
       height = 3)

# RPFs in function of mRNA ####
# no colorspace  
rpfvsmrna = abundNormLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = rna_ribofraction,
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("RPF Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks()

# saving
ggsave(filename = "plots/31_rpf_vs_mrna_tpwise.png",
       plot = rpfvsmrna,
       units = "in",
       width = 6,
       height = 5)

# Transcriptional Regulation in function of mRNA ####
# no colorspace
tcvsmrna = abundNormLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(protein_lysate / rna_total),
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Transcriptional Regulation") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
#  scale_y_log10(breaks = breaks,
#                minor_breaks = minor_breaks,
#                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b")

# saving
ggsave(filename = "plots/41_tc_vs_mrna_tpwise.png",
       plot = tcvsmrna,
       units = "in",
       width = 6,
       height = 5)

# slide tc from protein TPi+1 / mRNA TPi vs mRNA from TPi
# TP2 vs TP1
slide1 = abundNorm %>% 
  ggplot(aes(y = log2(mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP1),
             x = mean_abundance_rna_total_TP1)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Transcriptional Regulation (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
#  scale_y_log10(breaks = breaks,
#                minor_breaks = minor_breaks,
#                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b") +
  ggtitle("TC (TP2 / TP1) vs. mRNA TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = log2(mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP2),
             x = mean_abundance_rna_total_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Transcriptional Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
#  scale_y_log10(breaks = breaks,
#                minor_breaks = minor_breaks,
#                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b") +
  ggtitle("TC (TP3 / TP2) vs. mRNA TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = log2(mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP3),
             x = mean_abundance_rna_total_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Transcriptional Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
#  scale_y_log10(breaks = breaks,
#                minor_breaks = minor_breaks,
#                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b") +
  ggtitle("TC (TP4 / TP3) vs. mRNA TP3")

# arranging plots
tcvsmrnaslides = ggarrange(slide1, slide2, slide3,
                           ncol = 3, nrow = 1)

# saving
ggsave(filename = "plots/42_tc_vs_mrna_slides.png",
       plot = tcvsmrnaslides,
       units = "in",
       width = 10,
       height = 3)

# Translational Efficiency in function of mRNA ####
# no colorspace
tevsmrna = abundNormLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_ribofraction / rna_total),
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Translational Efficiency (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
#  scale_y_log10(breaks = breaks,
#                minor_breaks = minor_breaks,
#                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b")

# saving
ggsave(filename = "plots/51_te_vs_mrna_tpwise.png",
       plot = tevsmrna,
       units = "in",
       width = 6,
       height = 5)

# Translational Regulation (TLR) in function of mRNA ####
# no colorspace
tlrvsmrna = abundNormLongFuncat %>% 
  select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2((protein_lysate / rna_total) / (rna_ribofraction / rna_total)),
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  #  scale_y_log10(breaks = breaks,
  #                minor_breaks = minor_breaks,
  #                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b")

# saving
ggsave(filename = "plots/61_tlr_vs_mrna_tpwise.png",
       plot = tlrvsmrna,
       units = "in",
       width = 6,
       height = 5)

# slide tlr from TPi+1 vs mRNA from TPi
# TP2 vs TP1
slide1 = abundNorm %>% 
  ggplot(aes(y = log2((mean_abundance_protein_lysate_TP2 / mean_abundance_rna_total_TP1) / (mean_abundance_rna_ribofraction_TP1 / mean_abundance_rna_total_TP1)),
             x = mean_abundance_rna_total_TP1)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  #  scale_y_log10(breaks = breaks,
  #                minor_breaks = minor_breaks,
  #                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b") +
  ggtitle("TLR (TP2 / TP1) vs. mRNA TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = log2((mean_abundance_protein_lysate_TP3 / mean_abundance_rna_total_TP2) / (mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_total_TP2)),
             x = mean_abundance_rna_total_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  #  scale_y_log10(breaks = breaks,
  #                minor_breaks = minor_breaks,
  #                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b") +
  ggtitle("TLR (TP3 / TP2) vs. mRNA TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = log2((mean_abundance_protein_lysate_TP4 / mean_abundance_rna_total_TP3) / (mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_total_TP3)),
             x = mean_abundance_rna_total_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = TRUE)) +
  #  scale_y_log10(breaks = breaks,
  #                minor_breaks = minor_breaks,
  #                labels = function(x) format(x, scientific = TRUE)) +
  annotation_logticks(sides = "b") +
  ggtitle("TLR (TP4 / TP3) vs. mRNA TP3")

# arranging plots
tlrvsmrnaslides = ggarrange(slide1, slide2, slide3,
                            ncol = 3, nrow = 1)

# saving
ggsave(filename = "plots/62_tlr_vs_mrna_slides.png",
       plot = tlrvsmrnaslides,
       units = "in",
       width = 10,
       height = 3)

# trajectories of genes showing
# interesting patterns of translational regulation
melements = abundNormLongFuncat %>% 
  filter(arCOG == "Mobilome: prophages, transposons") %>% 
  select(locus_tag) %>% 
  unlist(use.names = F) %>% 
  unique()

# plotting all genes for manual inspection: abundance #####
allLocusTags = abundNormLongFuncat$locus_tag %>%
  sort() %>%
  unique()

# setting vars
i=1 # init
f=30 # final

# starting loop
while(f <= length(allLocusTags)){
  
  if(!dir.exists("plots/abundance")){
    dir.create("plots/abundance")
  }
  
  protSet = allLocusTags[i:f]
  
  # raw variables
  cols = c("protein_lysate" = "#E15759",
           "rna_ribofraction" = "#59A14F",
           "rna_total" = "#4E79A7")
  breaks = 10^(-10:10)
  minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  abundTrajectories = abundNormLongFuncat %>% 
    filter(locus_tag %in% protSet) %>% 
    filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
    ggplot(aes(x = timepoint,
               y = mean,
               color = libtype,
               group = libtype)) +
    geom_line() +
    geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
    facet_wrap(~ locus_tag) +
    scale_y_log10(breaks = breaks,
                  minor_breaks = minor_breaks,
                  labels = function(x) format(x, scientific = TRUE)) +
    ylab("Abundance") +
    xlab("Time point") +
    annotation_logticks(sides = "l") +
    scale_color_manual(name = "Lib. Type",
                       values = cols,
                       breaks = c("protein_lysate", "rna_ribofraction", "rna_total"),
                       labels = c("Protein", "RPFs", "mRNA"))
  
  # saving
  filename = paste0("plots/abundance/71_", i, "-", f, "abundTrajectories_tpwise.png")
  ggsave(filename = filename,
         plot = abundTrajectories,
         units = "in",
         width = 10,
         height = 8)
  
  if(!dir.exists("plots/ratios")){
    dir.create("plots/ratios")
  }
  
  # derived variables
  cols = c("rna_psiTE" = "#B07AA1",
           "rna_occupancy" = "#76B7B2",
           "rna_tlr" = "#F28E2B")
  ratioTrajectories = abundNormDerLong %>% 
    filter(locus_tag %in% protSet) %>% 
    filter(libtype != "protein_ribo") %>% 
    filter(libtype != "protein_lysate") %>% 
    filter(libtype != "rna_ribofraction") %>% 
    filter(libtype != "rna_total") %>% 
    ggplot(aes(x = timepoint,
               y = log2(mean),
               color = libtype,
               group = libtype)) +
    geom_line() +
    geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
    facet_wrap(~ locus_tag) +
    #  ylim(c(-12, 5)) +
    ylab("log2(Ratio)") +
    xlab("Time point") +
    scale_color_manual(name = "Ratio",
                       values = cols,
                       breaks = c("rna_occupancy", "rna_psiTE", "rna_tlr"),
                       labels = c("TE", "TC", "TLR"))
  
  # saving
  filename = paste0("plots/ratios/72_", i, "-", f, "ratioTrajectories_tpwise.png")
  ggsave(filename = filename,
         plot = ratioTrajectories,
         units = "in",
         width = 10,
         height = 8)
  
  # setting new vars
  i=f+1
  f=f+30
}

# saving individual plots as list items: abundance #####
trajectoryPlots = list()

for(i in allLocusTags){
  # raw variables
  cols = c("protein_lysate" = "#E15759",
           "rna_ribofraction" = "#59A14F",
           "rna_total" = "#4E79A7")
  breaks = 10^(-10:10)
  minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  trajectoryPlots[["abund"]][[i]] = abundNormLongFuncat %>% 
    filter(locus_tag %in% i) %>% 
    filter(libtype != "protein_ribo" & libtype != "rna_occupancy" & libtype != "rna_psiTE") %>% 
    ggplot(aes(x = timepoint,
               y = mean,
               color = libtype,
               group = libtype)) +
    geom_line() +
    geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
    facet_wrap(~ locus_tag) +
    scale_y_log10(breaks = breaks,
                  minor_breaks = minor_breaks,
                  labels = function(x) format(x, scientific = TRUE)) +
    ylab("Abundance") +
    xlab("Time point") +
    annotation_logticks(sides = "l") +
    scale_color_manual(name = "Lib. Type",
                       values = cols,
                       breaks = c("protein_lysate", "rna_ribofraction", "rna_total"),
                       labels = c("Protein", "RPFs", "mRNA"))
  
  # derived variables
  cols = c("rna_psiTE" = "#B07AA1",
           "rna_occupancy" = "#76B7B2",
           "rna_tlr" = "#F28E2B")
  trajectoryPlots[["ratios"]][[i]] = abundNormDerLong %>% 
    filter(locus_tag %in% i) %>% 
    filter(libtype != "protein_ribo") %>% 
    filter(libtype != "protein_lysate") %>% 
    filter(libtype != "rna_ribofraction") %>% 
    filter(libtype != "rna_total") %>% 
    ggplot(aes(x = timepoint,
               y = log2(mean),
               color = libtype,
               group = libtype)) +
    geom_line() +
    geom_linerange(aes(ymin = mean-se, ymax = mean+se)) +
    facet_wrap(~ locus_tag) +
    #  ylim(c(-12, 5)) +
    ylab("log2(Ratio)") +
    xlab("Time point") +
    scale_color_manual(name = "Ratio",
                       values = cols,
                       breaks = c("rna_occupancy", "rna_psiTE", "rna_tlr"),
                       labels = c("TE", "TC", "TLR"))
}


# plotting all genes for manual inspection: relative changes #####
# trajectories of protein, mRNA, and RPFs
allLocusTags = tc$locus_tag %>% 
  sort() %>% 
  unique()

# setting vars
i=1 # init
f=30 # final

# starting loop
while(f <= length(allLocusTags)){
  
  protSet = allLocusTags[i:f]
  
  if(!dir.exists("plots/lfc")){
    dir.create("plots/lfc")
  }
  
  cols = c("protein" = "#E15759",
           "RPF" = "#59A14F",
           "mRNA" = "#4E79A7")
  lfcTrajectories = tc %>%
    filter(tc$locus_tag %in% protSet) %>% 
    filter(timepoint != "TP32" & timepoint != "TP43") %>% 
    filter(libType != "RO" & libType != "beta" & libType != "TLR") %>% 
    ggplot(aes(x=timepoint, y=lfc, colour = libType, group = libType)) +
    geom_line() +
    geom_linerange(aes(ymin = lfc-lfcse, ymax = lfc+lfcse)) +
    facet_wrap(~locus_tag) +
    xlab("Time point") + ylab("log2(Fold Change)") +
    scale_color_manual(name = "Lib. Type",
                       values = cols,
                       breaks = c("protein", "RPF", "mRNA"),
                       labels = c("Protein", "RPFs", "mRNA"))
  
  # saving
  filename = paste0("plots/lfc/73_", i, "-", f, "lfcTrajectories_tpwise.png")
  ggsave(filename = filename,
         plot = lfcTrajectories,
         units = "in",
         width = 10,
         height = 8)
  
  if(!dir.exists("plots/lfcRatios")){
    dir.create("plots/lfcRatios")
  }
  
  # trajectories of ratios
  cols = c("beta" = "#B07AA1",
           "RO" = "#76B7B2",
           "TLR" = "#F28E2B")
  
  lfcRatioTrajectories = tc %>%
    filter(tc$locus_tag %in% protSet) %>% 
    filter(timepoint != "TP32" & timepoint != "TP43") %>% 
    filter(libType == "RO" | libType == "beta" | libType == "TLR") %>% 
    ggplot(aes(x=timepoint, y=lfc, colour = libType, group = libType)) +
    geom_line() +
    geom_linerange(aes(ymin = lfc-lfcse, ymax = lfc+lfcse)) +
    facet_wrap(~locus_tag) +
    xlab("Time point") + ylab("log2(Ratio)") +
    scale_color_manual(name = "Ratio",
                       values = cols,
                       breaks = c("RO", "beta", "TLR"),
                       labels = c("TE", "TC", "TLR"))
  
  # saving
  filename = paste0("plots/lfcRatios/74_", i, "-", f, "lfcRatioTrajectories_tpwise.png")
  ggsave(filename = filename,
         plot = lfcRatioTrajectories,
         units = "in",
         width = 10,
         height = 8)
  
  # setting new vars
  i=f+1
  f=f+30
}

# saving individual plots as list items: lfc #####
for(i in allLocusTags){
  # trajectories of foldchanges
  cols = c("protein" = "#E15759",
           "RPF" = "#59A14F",
           "mRNA" = "#4E79A7")
  trajectoryPlots[["lfc"]][[i]] = tc %>%
    filter(tc$locus_tag %in% i) %>% 
    filter(timepoint != "TP32" & timepoint != "TP43") %>% 
    filter(libType != "RO" & libType != "beta" & libType != "TLR") %>% 
    ggplot(aes(x=timepoint, y=lfc, colour = libType, group = libType)) +
    geom_line() +
    geom_linerange(aes(ymin = lfc-lfcse, ymax = lfc+lfcse)) +
    facet_wrap(~locus_tag) +
    xlab("Time point") + ylab("log2(Fold Change)") +
    scale_color_manual(name = "Lib. Type",
                       values = cols,
                       breaks = c("protein", "RPF", "mRNA"),
                       labels = c("Protein", "RPFs", "mRNA"))
  
  # trajectories of ratios
  cols = c("beta" = "#B07AA1",
           "RO" = "#76B7B2",
           "TLR" = "#F28E2B")
  
  trajectoryPlots[["lfcRatios"]][[i]] = tc %>%
    filter(tc$locus_tag %in% i) %>% 
    filter(timepoint != "TP32" & timepoint != "TP43") %>% 
    filter(libType == "RO" | libType == "beta" | libType == "TLR") %>% 
    ggplot(aes(x=timepoint, y=lfc, colour = libType, group = libType)) +
    geom_line() +
    geom_linerange(aes(ymin = lfc-lfcse, ymax = lfc+lfcse)) +
    facet_wrap(~locus_tag) +
    xlab("Time point") + ylab("log2(Ratio)") +
    scale_color_manual(name = "Ratio",
                       values = cols,
                       breaks = c("RO", "beta", "TLR"),
                       labels = c("TE", "TC", "TLR"))
}

# a list of manually inspected genes: proteins explained by mRNA #####
# here, protein levels are explained by mrna levels
all = c("VNG_0614G",
        "VNG_0013C",
        "VNG_0084G")

ggsave(filename = "plots/71_abundTrajectories_tc_tpwise.png",
       plot = abundTrajectories,
       units = "in",
       width = 10,
       height = 3)

# a list of manually inspected genes: proteins explained by occupancy #####
# here, protein levels are explained by mrna levels and occupancy
all = c("VNG_0063G",
        "VNG_0115G",
        "VNG_1204G")

ggsave(filename = "plots/71_abundTrajectories_te_tpwise.png",
       plot = abundTrajectories,
       units = "in",
       width = 10,
       height = 3)
