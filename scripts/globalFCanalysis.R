# 20210308 alorenzetti

# description ####
# this script will take the non normalized
# abundance values and compute time point-wise
# fold changes
# then the genes are going to be grouped
# according to the direction of change

# vars ####
# setting log2 FC minimum difference
# to change status
fcthr = 1.5
lfcthr = log2(1.5)

# if user wants NCBI locus_tags
# as primary, change this var to
# y
# locus_tag NCBI
ncbiLocusTag = "y"

# saving NCBI locus_tag dictionary
ncbiDict = nrtx %>% 
  separate_rows(locus_tag, sep = ",") %>% 
  dplyr::filter(str_detect(locus_tag, "VNG_RS")) %>% 
  group_by(representative) %>% 
  summarise(locus_tag = paste0(locus_tag, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(locus_tag = str_replace(locus_tag, ",.*$", ""))

# starting processing ####
# creating a new object to compute fold changes
# and changing
fcs = abund %>% 
  mutate(mean_lfc_rna_total_TP2vsTP1 = log2(mean_abundance_rna_total_TP2 / mean_abundance_rna_total_TP1),
         mean_lfc_rna_total_TP3vsTP2 = log2(mean_abundance_rna_total_TP3 / mean_abundance_rna_total_TP2),
         mean_lfc_rna_total_TP4vsTP3 = log2(mean_abundance_rna_total_TP4 / mean_abundance_rna_total_TP3),
         
         mean_lfc_rna_ribofraction_TP2vsTP1 = log2(mean_abundance_rna_ribofraction_TP2 / mean_abundance_rna_ribofraction_TP1),
         mean_lfc_rna_ribofraction_TP3vsTP2 = log2(mean_abundance_rna_ribofraction_TP3 / mean_abundance_rna_ribofraction_TP2),
         mean_lfc_rna_ribofraction_TP4vsTP3 = log2(mean_abundance_rna_ribofraction_TP4 / mean_abundance_rna_ribofraction_TP3),
         
         mean_lfc_protein_lysate_TP2vsTP1 = log2(mean_abundance_protein_lysate_TP2 / mean_abundance_protein_lysate_TP1),
         mean_lfc_protein_lysate_TP3vsTP2 = log2(mean_abundance_protein_lysate_TP3 / mean_abundance_protein_lysate_TP2),
         mean_lfc_protein_lysate_TP4vsTP3 = log2(mean_abundance_protein_lysate_TP4 / mean_abundance_protein_lysate_TP3)) %>% 
  select(c(locus_tag, contains("lfc")))

# replacing locus_tags if requested by the user
if(ncbiLocusTag == "y"){
  fcs = fcs %>% 
    left_join(x = ., y = ncbiDict,
              by = c("locus_tag" = "representative")) %>% 
    dplyr::mutate(locus_tag = locus_tag.y) %>% 
    dplyr::select(-locus_tag.y) %>% 
    filter(!is.na(locus_tag))
}

# creating new columns to store change status
fcs = fcs %>% 
  mutate(change_rna_total_TP2vsTP1 = case_when(mean_lfc_rna_total_TP2vsTP1 > lfcthr ~ "Up",
                                               mean_lfc_rna_total_TP2vsTP1 < -lfcthr ~ "Down",
                                               abs(mean_lfc_rna_total_TP2vsTP1) < lfcthr ~ "Flat",
                                               TRUE ~ "?"),
         change_rna_total_TP3vsTP2 = case_when(mean_lfc_rna_total_TP3vsTP2 > lfcthr ~ "Up",
                                               mean_lfc_rna_total_TP3vsTP2 < -lfcthr ~ "Down",
                                               abs(mean_lfc_rna_total_TP3vsTP2) < lfcthr ~ "Flat",
                                               TRUE ~ "?"),
         change_rna_total_TP4vsTP3 = case_when(mean_lfc_rna_total_TP4vsTP3 > lfcthr ~ "Up",
                                               mean_lfc_rna_total_TP4vsTP3 < -lfcthr ~ "Down",
                                               abs(mean_lfc_rna_total_TP4vsTP3) < lfcthr ~ "Flat",
                                               TRUE ~ "?")) %>% 
  
  mutate(change_rna_ribofraction_TP2vsTP1 = case_when(mean_lfc_rna_ribofraction_TP2vsTP1 > lfcthr ~ "Up",
                                                      mean_lfc_rna_ribofraction_TP2vsTP1 < -lfcthr ~ "Down",
                                                      abs(mean_lfc_rna_ribofraction_TP2vsTP1) < lfcthr ~ "Flat",
                                                      TRUE ~ "?"),
         change_rna_ribofraction_TP3vsTP2 = case_when(mean_lfc_rna_ribofraction_TP3vsTP2 > lfcthr ~ "Up",
                                                      mean_lfc_rna_ribofraction_TP3vsTP2 < -lfcthr ~ "Down",
                                                      abs(mean_lfc_rna_ribofraction_TP3vsTP2) < lfcthr ~ "Flat",
                                                      TRUE ~ "?"),
         change_rna_ribofraction_TP4vsTP3 = case_when(mean_lfc_rna_ribofraction_TP4vsTP3 > lfcthr ~ "Up",
                                                      mean_lfc_rna_ribofraction_TP4vsTP3 < -lfcthr ~ "Down",
                                                      abs(mean_lfc_rna_ribofraction_TP4vsTP3) < lfcthr ~ "Flat",
                                                      TRUE ~ "?")) %>% 
  
  mutate(change_protein_lysate_TP2vsTP1 = case_when(mean_lfc_protein_lysate_TP2vsTP1 > lfcthr ~ "Up",
                                                    mean_lfc_protein_lysate_TP2vsTP1 < -lfcthr ~ "Down",
                                                    abs(mean_lfc_protein_lysate_TP2vsTP1) < lfcthr ~ "Flat",
                                                    TRUE ~ "?"),
         change_protein_lysate_TP3vsTP2 = case_when(mean_lfc_protein_lysate_TP3vsTP2 > lfcthr ~ "Up",
                                                    mean_lfc_protein_lysate_TP3vsTP2 < -lfcthr ~ "Down",
                                                    abs(mean_lfc_protein_lysate_TP3vsTP2) < lfcthr ~ "Flat",
                                                    TRUE ~ "?"),
         change_protein_lysate_TP4vsTP3 = case_when(mean_lfc_protein_lysate_TP4vsTP3 > lfcthr ~ "Up",
                                                    mean_lfc_protein_lysate_TP4vsTP3 < -lfcthr ~ "Down",
                                                    abs(mean_lfc_protein_lysate_TP4vsTP3) < lfcthr ~ "Flat",
                                                    TRUE ~ "?")) %>% 
  
  # creating a code for the time course changes
  mutate(timecourse_code_rna_total = paste(change_rna_total_TP2vsTP1,
                                           change_rna_total_TP3vsTP2,
                                           change_rna_total_TP4vsTP3,
                                           sep = "_"),
         timecourse_code_rna_ribofraction = paste(change_rna_ribofraction_TP2vsTP1,
                                                  change_rna_ribofraction_TP3vsTP2,
                                                  change_rna_ribofraction_TP4vsTP3,
                                                  sep = "_"),
         timecourse_code_protein_lysate = paste(change_protein_lysate_TP2vsTP1,
                                                change_protein_lysate_TP3vsTP2,
                                                change_protein_lysate_TP4vsTP3,
                                                sep = "_"))

# organizing data structure to get insights
# about the genes following each one of the 
# change codes
tcCodes = list()

# everything
tcCodes[["all"]] = fcs %>%
  select(c(locus_tag,
           starts_with("timecourse_code")))

# subsetting and grouping rna, ribo, and protein
# rna
tcCodes[["rna_total"]] = tcCodes[["all"]] %>%
  select(c(locus_tag,
           contains("rna_total"))) %>% 
  group_by(timecourse_code_rna_total) %>% 
  summarise(entries_rna_total = paste0(locus_tag, collapse = ","))

# ribo 
tcCodes[["rna_ribofraction"]] = tcCodes[["all"]] %>%
  select(c(locus_tag,
           contains("rna_ribofraction"))) %>% 
  group_by(timecourse_code_rna_ribofraction) %>% 
  summarise(entries_rna_ribofraction = paste0(locus_tag, collapse = ","))

# protein
tcCodes[["protein_lysate"]] = tcCodes[["all"]] %>%
  select(c(locus_tag,
           contains("protein_lysate"))) %>% 
  group_by(timecourse_code_protein_lysate) %>% 
  summarise(entries_protein_lysate = paste0(locus_tag, collapse = ","))

# unifying datasets 
tcCodes[["unified_written"]] = full_join(x = tcCodes[["rna_total"]],
                                         y = tcCodes[["rna_ribofraction"]],
                                         by = c("timecourse_code_rna_total" = "timecourse_code_rna_ribofraction"))

tcCodes[["unified_written"]] = full_join(x = tcCodes[["unified_written"]],
                                         y = tcCodes[["protein_lysate"]],
                                         by = c("timecourse_code_rna_total" = "timecourse_code_protein_lysate")) %>% 
  dplyr::rename(timecourse_code = "timecourse_code_rna_total")

# ordering data set
tcCodes[["unified_written"]] = tcCodes[["unified_written"]][mixedorder(tcCodes[["unified_written"]]$timecourse_code,
                                                                       decreasing = T),]
# creating cols to store intersections
# and then finding them
tcCodes$unified_written$intersection_total_prot = NA_character_
tcCodes$unified_written$intersection_total_ribo = NA_character_
tcCodes$unified_written$intersection_prot_ribo = NA_character_
tcCodes$unified_written$intersection_all = NA_character_

for(i in 1:dim(tcCodes$unified_written)[1]){
  
  rna = tcCodes$unified_written$entries_rna_total[i] %>%
    str_split(pattern = ",", simplify = F) %>%
    unlist()
  
  ribo = tcCodes$unified_written$entries_rna_ribofraction[i] %>%
    str_split(pattern = ",", simplify = F) %>%
    unlist()
  
  prot = tcCodes$unified_written$entries_protein_lysate[i] %>%
    str_split(pattern = ",", simplify = F) %>%
    unlist()
  
  rnaprot = base::intersect(rna, prot) ; if(identical(rnaprot, character(0))){rnaprot = NA_character_}
  rnaribo = base::intersect(rna, ribo) ; if(identical(rnaribo, character(0))){rnaribo = NA_character_}
  protribo = base::intersect(prot, ribo) ; if(identical(protribo, character(0))){protribo = NA_character_}
  all = base::intersect(base::intersect(rna, prot), ribo) ; if(identical(all, character(0))){all = NA_character_}
  
  tcCodes$unified_written$intersection_total_prot[i] = paste0(rnaprot, collapse = ",")
  tcCodes$unified_written$intersection_total_ribo[i] = paste0(rnaribo, collapse = ",")
  tcCodes$unified_written$intersection_prot_ribo[i] = paste0(protribo, collapse = ",")
  tcCodes$unified_written$intersection_all[i] = paste0(all, collapse = ",")
}

# adjusting NAs of intersection cols
# and headers
tcCodes[["unified_written"]] = tcCodes[["unified_written"]] %>% 
  mutate(across(.cols = starts_with("intersection"),
                .fns = ~ case_when(.x  == "NA" ~ NA_character_,
                                   TRUE ~ as.character(.x))))

# renaming cols
colnames(tcCodes[["unified_written"]]) = c("timecourse_code",
                                           "Entries mRNA",
                                           "Entries RiboSeq",
                                           "Entries Protein",
                                           "Intersection mRNA Protein",
                                           "Intersection mRNA RiboSeq",
                                           "Intersection Protein RiboSeq",
                                           "Intersection mRNA Protein RiboSeq")

# creating a count matrix based in the above matrix
tcCodes[["unified_counts"]] = tcCodes[["unified_written"]] %>% 
  mutate(across(.cols = -timecourse_code,
                .fns = ~ str_count(string = .x, pattern = "VNG_")))

tcCodes[["webObj"]] = map2_df(tcCodes[["unified_counts"]],
                              tcCodes[["unified_written"]],
                              ~ paste(.x, .y)) %>% 
  mutate(timecourse_code = str_replace(timecourse_code, " .*$", "")) %>% 
  mutate(across(.cols = -timecourse_code,
                .fns = ~ case_when(str_detect(.x, "NA NA") ~ NA_character_,
                                   TRUE ~ str_replace(string = .x,
                                                      pattern = '^(.*) (.*)$',
                                                      replacement = '<a href="#" onclick="alert(\'\\2\');">\\1</a>'))))  %>% 
  mutate(across(.cols = -timecourse_code,
                .fns = ~ case_when(is.na(.x) ~ NA_character_,
                                   TRUE ~ str_replace_all(string = .x,
                                                          pattern = ',',
                                                          replacement = ' '))))


tcCodes[["unified_written"]] %>% 
  write_tsv(x = .,
            file = "results/20210311_globalFCanalysis_entries_table.tsv")

tcCodes[["unified_counts"]] %>% 
  write_tsv(x = .,
            file = "results/20210311_globalFCanalysis_counts_table.tsv")

# generating and saving datatable obj
tcCodes[["webObj"]] %>% 
  datatable(data = ., escape = F) %>% 
  DT::saveWidget(widget = .,
                 file = "20210311_globalFCanalysis_int_table.html",
                 selfcontained = T)
