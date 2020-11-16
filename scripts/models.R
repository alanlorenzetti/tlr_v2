# alorenze 20200525

# description ####
# this script will generate
# linear models to explain
# proteins in function of mRNA and/or RPFs
# it will also plot several charts

# loading libs #####
source("scripts/loadingLibs.R")

# setting thresholds
log2fcthreshold = 0.75
borderlinezero = 0.25

# defining tab10 color scheme
# ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`
tab10col = c(
  "Q1"="#E15759", # red
  "Q2"="#F28E2B", # orange
  "Q3"="#4E79A7", # blue
  "Q4"="#59A14F", # green
  "BL12"="#FF9DA7", # pink
  "BL14"="#B07AA1", # purple
  "BL23"="#76B7B2", # light teal 
  "BL34"="#9C755F", # brown
  "core"="#EDC948", # yellow
  "center"="grey80") # grey

# mod to allow black lines using fill aes
tab10colmod = tab10col
tab10colmod["Q1"] = "black"
tab10colmod["Q3"] = "black"

# setting ggplot blank theme
theme_set(theme_bw())

# unifying datasets
alldfs = list()
quantileNorm = list()
for(i in 2:4 %>% as.character()){
  varName[1] = paste0("sigProt")
  varName[2] = paste0("borderlineZeroProt")
  varName[3] = paste0("sigtotrna")
  varName[4] = paste0("sigribo")
  varName[5] = paste0("borderlineZerototrna")
  varName[6] = paste0("borderlineZeroribo")
  
  alldfs[[paste0("TP",i)]] = inner_join(tpivstp1[[paste0("TP",i)]],
                                        unifiedFin[[paste0("TP",i)]],
                                        by = c("locus_tag" = "locus_tag")) %>% 
    dplyr::select(locus_tag,
           matches("lfc"),
           matches("sig"),
           matches("border")) %>%
    rename("mRNA" = mean_lfc_rna_total,
           "RPF" = mean_lfc_rna_ribo,
           "protein" = mean_lfc_protein_lysate,
           "lfcse_mRNA" = se_lfc_rna_total,
           "lfcse_RPF" = se_lfc_rna_ribo)
    
  # I will apply quantile normalization to mRNA, RPF and protein
  # and then I will replace the original values by the normalized ones
  # the normalized values correlate quite well with the original ones
  # so I don't see a reason to be concerned about adjusting significance
  # normalizing lfc values
  quantileNorm[[paste0("TP",i)]] = alldfs[[paste0("TP",i)]] %>%
    dplyr::select(mRNA, RPF, protein) %>%
    as.matrix() %>% 
    normalize.quantiles() %>% 
    as_tibble()
  colnames(quantileNorm[[paste0("TP",i)]]) = c("mRNA", "RPF", "protein")
  
  alldfs[[paste0("TP",i)]][,c("mRNA","RPF","protein")] = quantileNorm[[paste0("TP",i)]][,c("mRNA","RPF","protein")]
  
  # normalizing standard errors
  quantileNorm[[paste0("TP",i)]] = alldfs[[paste0("TP",i)]] %>%
    dplyr::select(lfcse_mRNA, lfcse_RPF) %>%
    as.matrix() %>% 
    normalize.quantiles() %>% 
    as_tibble()
  colnames(quantileNorm[[paste0("TP",i)]]) = c("lfcse_mRNA", "lfcse_RPF")
  
  alldfs[[paste0("TP",i)]][,c("lfcse_mRNA","lfcse_RPF")] = quantileNorm[[paste0("TP",i)]][,c("lfcse_mRNA","lfcse_RPF")]
  
  # finding significant and borderline status again
  # after quantile normalization
  for(type in c("totrna", "ribo", "Prot")){
    var1 = paste0("sig",type)
    var2 = paste0("borderlineZero",type)
    
    if(type == "totrna"){var3 = "mRNA"}
    if(type == "ribo"){var3 = "RPF"}
    if(type == "Prot"){var3 = "protein"}
    
    alldfs[[paste0("TP",i)]] = alldfs[[paste0("TP",i)]] %>% 
      mutate(!!var1 := case_when(abs(get(var3)) >= log2fcthreshold & get(var1) == "yes" ~ "yes",
                                 TRUE ~ "no")) %>% 
      mutate(!!var2 := case_when(abs(get(var3)) < borderlinezero ~ "yes",                                     TRUE ~ "no"))
  }
  
  # creating new categorical variables
  alldfs[[paste0("TP",i)]] = alldfs[[paste0("TP",i)]] %>% 
    mutate(`RPF-mRNA-quad` = case_when(mRNA >= log2fcthreshold & RPF >= log2fcthreshold ~ "Q1",
                                       mRNA <= -log2fcthreshold & RPF >= log2fcthreshold ~ "Q2",
                                       mRNA <= -log2fcthreshold & RPF <= -log2fcthreshold ~ "Q3",
                                       mRNA >= log2fcthreshold & RPF <= -log2fcthreshold ~ "Q4",
                                       mRNA >= borderlinezero & abs(RPF) < borderlinezero ~ "BL14",
                                       mRNA <= -borderlinezero & abs(RPF) < borderlinezero ~ "BL23",
                                       RPF >= borderlinezero & abs(mRNA) < borderlinezero ~ "BL12",
                                       RPF <= -borderlinezero & abs(mRNA) < borderlinezero ~ "BL34",
                                       TRUE ~ "center")) %>% 
    
    mutate(`protein-mRNA-quad` = case_when(mRNA >= log2fcthreshold & protein >= log2fcthreshold ~ "Q1",
                                           mRNA <= -log2fcthreshold & protein >= log2fcthreshold ~ "Q2",
                                           mRNA <= -log2fcthreshold & protein <= -log2fcthreshold ~ "Q3",
                                           mRNA >= log2fcthreshold & protein <= -log2fcthreshold ~ "Q4",
                                           mRNA >= borderlinezero & abs(protein) < borderlinezero ~ "BL14",
                                           mRNA <= -borderlinezero & abs(protein) < borderlinezero ~ "BL23",
                                           protein >= borderlinezero & abs(mRNA) < borderlinezero ~ "BL12",
                                           protein <= -borderlinezero & abs(mRNA) < borderlinezero ~ "BL34",
                                           TRUE ~ "center")) %>%
    
    mutate(`protein-RPF-quad` = case_when(RPF >= log2fcthreshold & protein >= log2fcthreshold ~ "Q1",
                                          RPF <= -log2fcthreshold & protein >= log2fcthreshold ~ "Q2",
                                          RPF <= -log2fcthreshold & protein <= -log2fcthreshold ~ "Q3",
                                          RPF >= log2fcthreshold & protein <= -log2fcthreshold  ~ "Q4",
                                          RPF >= borderlinezero & abs(protein) < borderlinezero ~ "BL14",
                                          RPF <= -borderlinezero & abs(protein) < borderlinezero ~ "BL23",
                                          protein >= borderlinezero & abs(RPF) < borderlinezero ~ "BL12",
                                          protein <= -borderlinezero & abs(RPF) < borderlinezero ~ "BL34",
                                          TRUE ~ "center")) %>% 
    
    mutate(`RPF-mRNA-quad` = case_when(get(varName[4]) == "yes" & get(varName[3]) == "yes" ~ as.character(`RPF-mRNA-quad`),
                                       get(varName[4]) == "yes" & get(varName[3]) == "no" ~ as.character(`RPF-mRNA-quad`),
                                       get(varName[4]) == "no" & get(varName[3]) == "yes" ~ as.character(`RPF-mRNA-quad`),
                                       abs(RPF) < borderlinezero & abs(mRNA) < borderlinezero ~ "core",
                                       TRUE ~ "center"),
           `protein-mRNA-quad` = case_when(get(varName[1]) == "yes" & get(varName[3]) == "yes" ~ as.character(`protein-mRNA-quad`),
                                           get(varName[1]) == "yes" & get(varName[3]) == "no" ~ as.character(`protein-mRNA-quad`),
                                           get(varName[1]) == "no" & get(varName[3]) == "yes" ~ as.character(`protein-mRNA-quad`),
                                           abs(protein) < borderlinezero & abs(mRNA) < borderlinezero ~ "core",
                                           TRUE ~ "center"),
           `protein-RPF-quad` = case_when(get(varName[1]) == "yes" & get(varName[4]) == "yes" ~ as.character(`protein-RPF-quad`),
                                          get(varName[1]) == "yes" & get(varName[4]) == "no" ~ as.character(`protein-RPF-quad`),
                                          get(varName[1]) == "no" & get(varName[4]) == "yes" ~ as.character(`protein-RPF-quad`),
                                          abs(protein) < borderlinezero & abs(RPF) < borderlinezero ~ "core",
                                          TRUE ~ "center")) %>% 
    
    mutate(RPF_mRNA_sig = case_when(get(varName[4]) == "yes" & get(varName[3]) == "yes" ~ "bold",
                                    get(varName[4]) == "yes" & get(varName[3]) == "no" ~ "mid",
                                    get(varName[4]) == "no" & get(varName[3]) == "yes" ~ "mid",
                                    TRUE ~ "faint"),
           protein_mRNA_sig = case_when(get(varName[1]) == "yes" & get(varName[3]) == "yes" ~ "bold",
                                        get(varName[1]) == "yes" & get(varName[3]) == "no" ~ "mid",
                                        get(varName[1]) == "no" & get(varName[3]) == "yes" ~ "mid",
                                        TRUE ~ "faint"),
           protein_RPF_sig = case_when(get(varName[1]) == "yes" & get(varName[4]) == "yes" ~ "bold",
                                       get(varName[1]) == "yes" & get(varName[4]) == "no" ~ "mid",
                                       get(varName[1]) == "no" & get(varName[4]) == "yes" ~ "mid",
                                       TRUE ~ "faint")) %>% 
    
    mutate(timePoint = i %>% as.numeric())
}

# functions to find axes lims
findXLims = function(df){
  mrnamax = df[,"mRNA"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  rpfmax = df[,"RPF"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  xmax = max(mrnamax,rpfmax)
  xmin = -xmax
  Xlims = c(xmin,xmax)
  return(Xlims)
}
findYLims = function(df){
  rpfmax = df[,"RPF"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  protmax = df[,"protein"] %>% abs() %>% max() %>% RoundTo(x = ., multiple = 2, FUN = ceiling)
  ymax = max(rpfmax,protmax)
  ymin = -ymax
  Ylims = c(ymin,ymax)
  return(Ylims)
}

# getting models to explain relationship
# between y and x
modelGeneral = function(df, type){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"}
  md = lm(formula=paste0(y,"~",x),data=df)
  return(md)
}

# running modelling
models = list()
types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")
for(i in names(alldfs)){
  for(k in types){
    models[[i]][[k]] = modelGeneral(alldfs[[i]], k) %>% summary()
  }
}

# plotting general trends
plotGeneral = function(df, type, xLim, yLim){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x)}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x)}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x)}
  
  model = lm(formula = form, data = df)
  cortext = paste0("*R<sup>2</sup>* = ", summary(model)$r.squared %>% round(digits = 3))
  p = pf(summary(model)$fstatistic[1],
         summary(model)$fstatistic[2],
         summary(model)$fstatistic[3],
         lower.tail = FALSE) %>% 
    formatC(format = "e", digits = 2)
  pvaltext = paste0("*p* = ", p)
  
  grob = grobTree(richtext_grob(text = paste0(cortext, "; ", pvaltext),
                           gp=gpar(fontsize=8),
                           x=0.075, y=0.95, hjust=0))
  
  alp = 0.05
  sz = 2
  
  plot = ggplot(data = df) +
    geom_point(aes(x=get(x), y=get(y)), alpha = alp) +
    xlim(xLim) + ylim(yLim) + xlab(x) + ylab(y) +
#    geom_vline(xintercept = 0, size = 0.3) +
#    geom_hline(yintercept = 0, size = 0.3) +
    geom_point(inherit.aes = F,
                data = df,
                aes(x=get(x), y=get(y)),
                alpha = alp,
                show.legend = F) +
    geom_smooth(inherit.aes = F,
                data = df,
                aes(x=get(x), y=get(y)),
                method="lm",
                alpha = alp,
                color = "blue",
                show.legend = F,
                size = 0.5) +
    annotation_custom(grob = grob)
  
  return(plot)
}

# generating general plots
genplots = list()
for(i in names(alldfs)){
  xLim = findXLims(alldfs[[i]])
  yLim = findYLims(alldfs[[i]])
  for(k in types){
    genplots[[i]][[k]] = plotGeneral(alldfs[[i]], k, xLim, yLim)}
}

# arranging plots (full)
ggarrange(plotlist = c(genplots$TP2,
                       genplots$TP3,
                       genplots$TP4),
          nrow = 3, ncol = 3)

# plot functions with sig dataset
plotLM = function(df, type, xLim, yLim){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}

  if(y == "RPF"){var1 = "sigribo"}
  if(y == "protein"){var1 = "sigProt"}
  if(x == "mRNA"){var2 = "sigtotrna"}
  if(x == "RPF"){var2 = "sigribo"}
  
  varName1 = df %>%
    dplyr::select(starts_with(var1)) %>% names() %>% as.symbol()
  varName2 = df %>%
    dplyr::select(starts_with(var2)) %>% names() %>% as.symbol()
  colsymb = color %>% as.symbol()
  dfSig = df %>% filter((!!varName1 == "yes" & !!varName2 == "yes") & (!!colsymb == "Q1" | !!colsymb == "Q3"))
#  dfNonSig =  df %>% filter(!!varName2 == "no")
  
#  alp = 0.25
  sz = 1.5
  
  var = color %>% as.symbol()
  model = list() ; cortext = list()
  p = list() ; pvaltext = list()
  grob = list() ; grobY = 0.92
  
  quadsItSummy = dfSig %>% dplyr::select(!!var) %>% table()
  quadsIt = names(quadsItSummy[quadsItSummy > 1])
  quadsIt = quadsIt[quadsIt %in% paste0("Q", c(1,3))]
  
  for(i in quadsIt){
    model[[i]] = lm(formula = form, data = dfSig %>% filter(!!var == i))
    cortext[[i]] = paste0("*R<sup>2</sup>* = ", summary(model[[i]])$r.squared %>% round(digits = 3))
    p[[i]] = pf(summary(model[[i]])$fstatistic[1],
                summary(model[[i]])$fstatistic[2],
                summary(model[[i]])$fstatistic[3],
                lower.tail = FALSE) %>% 
      formatC(format = "e", digits = 2)
    pvaltext[[i]] = paste0("*p* = ", p[[i]])
    if(i == "Q1"){col = tab10col["Q1"]}
#    if(i == "Q2"){col = tab10col["Q2"]}
    if(i == "Q3"){col = tab10col["Q3"]}
#    if(i == "Q4"){col = tab10col["Q4"]}
    grob[[i]] = grobTree(richtext_grob(text = paste0(cortext[[i]], "; ", pvaltext[[i]]),
                                  gp=gpar(fontsize=8, col=col),
                                  x=0.075, y=grobY, hjust=0))
    grobY = grobY - 0.05
  }
  
  plot = ggplot(data = df, aes(x=get(x), y=get(y), colour=get(color))) +
    geom_point(aes(alpha=get(alpha)), size = sz, show.legend = F) +
    xlim(xLim) + ylim(yLim) + xlab(x) + ylab(y) +
    # geom_vline(xintercept = 0, size = 0.3) +
    # geom_hline(yintercept = 0, size = 0.3) +
    # geom_point(inherit.aes = F,
    #             data = dfSig,
    #             aes(x=get(x), y=get(y), color=get(color)),
    #             alpha = alp, size = sz,
    #             show.legend = F) +
    geom_smooth(inherit.aes = F,
                data = dfSig,
                aes(x=get(x), y=get(y), fill=get(color)),
                colour="black",
                size=0.5,
                method="lm",
                show.legend = F) +
    scale_colour_manual(values = tab10col) +
    scale_fill_manual(values = tab10col) +
    scale_alpha_manual(values = c("bold"=0.8, "mid"= 0.3, "faint"=0.1))
#    scale_shape_manual(values = c("bold"=4, "mid"=5, "faint"=1))
  
  for(i in quadsIt){
    plot = plot +
      annotation_custom(grob = grob[[i]])
  }
  
  return(plot)
}

# types of plot we should get
types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")

######## linear models and plots for full dataset
lmplots = list()
for(i in names(alldfs)){
  xLim = findXLims(alldfs[[i]])
  yLim = findYLims(alldfs[[i]])
  for(k in types){
    lmplots[[i]][[k]] = plotLM(alldfs[[i]], k, xLim, yLim)
  }
}

# arranging plots (full panel)
ggarrange(plotlist = c(lmplots$TP2,
                       lmplots$TP3,
                       lmplots$TP4),
          nrow = 3, ncol = 3)
