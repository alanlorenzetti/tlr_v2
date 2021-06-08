# 20201009
# request was made on 20201009

# description #####
# Nitin has requested a couple of charts to show the correlation between 
# transcripts, RPF, and protein levels. Ideally, we should be able to
# identify cases of transcriptional regulation, translational efficiency
# and translational regulation

# plotting ####
breaks = 10^(-10:10)
minor_breaks = rep(1:9, 21)*(10^rep(-10:10, each=9))

# protein abundance in function of mRNA ####
# no colorspace  
protvsmrna = abundNormLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  # ylab("Abundância de Proteínas") +
  # xlab("Abundância de mRNAs") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. mRNA TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_total_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. mRNA TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_total_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
proteinvsrpf = abundNormLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = protein_lysate,
             x = rna_ribofraction)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Protein Abundance") +
  xlab("RPF Abundance") +
  # ylab("Abundância de Proteínas") +
  # xlab("Abundância de RPFs") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("RPF Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP2 vs. RPF TP1")

# TP3 vs TP2
slide2 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP3,
             x = mean_abundance_rna_ribofraction_TP2)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("RPF Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks() +
  ggtitle("Protein TP3 vs. RPF TP2")

# TP4 vs TP3
slide3 = abundNorm %>% 
  ggplot(aes(y = mean_abundance_protein_lysate_TP4,
             x = mean_abundance_rna_ribofraction_TP3)) +
  geom_point(alpha = 0.25) + 
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Protein Abundance") +
  xlab("RPF Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
rpfvsmrna = abundNormLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = rna_ribofraction,
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  # ylab("RPF Abundance") +
  # xlab("mRNA Abundance") +
  ylab("Abundância de RPFs") +
  xlab("Abundância de mRNAs") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  scale_y_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
  annotation_logticks()

# saving
ggsave(filename = "plots/31_rpf_vs_mrna_tpwise.png",
       plot = rpfvsmrna,
       units = "in",
       width = 6,
       height = 5)

# Transcriptional Regulation in function of mRNA ####
# no colorspace
tcvsmrna = abundNormLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(protein_lysate / rna_total),
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Transcriptional Regulation") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Transcriptional Regulation (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Transcriptional Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Transcriptional Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
tevsmrna = abundNormLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2(rna_ribofraction / rna_total),
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Translational Efficiency (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
tlrvsmrna = abundNormLong %>% 
  filter(timepoint != "TP0") %>% 
  dplyr::select(-se) %>% 
  pivot_wider(names_from = libtype,
              values_from = mean) %>% 
  ggplot(aes(y = log2((protein_lysate / rna_total) / (rna_ribofraction / rna_total)),
             x = rna_total)) +
  #    geom_density2d() +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  facet_wrap(~ timepoint) +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance") +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  geom_smooth(method = "lm",
              color = tab10$value[1]) +
  stat_cor(method = "pearson") +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
  stat_cor(method = "pearson",
           color = tab10$value[1]) +
  ylab("Translational Regulation (log2)") +
  xlab("mRNA Abundance")  +
  scale_x_log10(breaks = breaks,
                minor_breaks = minor_breaks,
                labels = function(x) format(x, scientific = F)) +
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
