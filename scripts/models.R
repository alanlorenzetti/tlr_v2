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
  varName[1] = paste0("sigProtOnTP",i)
  varName[2] = paste0("borderlineZeroProtOnTP",i)
  varName[3] = paste0("sigtotrnaOnTP",i)
  varName[4] = paste0("sigriboOnTP",i)
  varName[5] = paste0("borderlineZerototrnaOnTP",i)
  varName[6] = paste0("borderlineZeroriboOnTP",i)
  
  alldfs[[paste0("TP",i)]] = inner_join(tpivstp1[[paste0("TP",i)]],
                                        unifiedFin[[paste0("TP",i)]],
                                        by = c("locus_tag" = "locus_tag")) %>% 
    select(locus_tag,
           log2FoldChange.x,
           lfcSE.x,
           log2FoldChange.y,
           lfcSE.y,
           lfc,
           lfcse,
           !!varName[1],
           !!varName[2],
           !!varName[3],
           !!varName[4],
           !!varName[5],
           !!varName[6]) %>%
    rename("mRNA" = log2FoldChange.x,
           "RPF" = log2FoldChange.y,
           "protein" = lfc,
           "lfcse_mRNA" = lfcSE.x,
           "lfcse_RPF" = lfcSE.y,
           "lfcse_protein" = lfcse)
    
  # I will apply quantile normalization to mRNA, RPF and protein
  # and then I will replace the original values by the normalized ones
  # the normalized values correlate quite well with the original ones
  # so I don't see a reason to be concerned about adjusting significance
  quantileNorm[[paste0("TP",i)]] = alldfs[[paste0("TP",i)]] %>%
    select(mRNA, RPF, protein) %>%
    as.matrix() %>% 
    normalize.quantiles() %>% 
    as_tibble()
  colnames(quantileNorm[[paste0("TP",i)]]) = c("mRNA", "RPF", "protein")
  
  alldfs[[paste0("TP",i)]][,c("mRNA","RPF","protein")] = quantileNorm[[paste0("TP",i)]][,c("mRNA","RPF","protein")]
  
  # finding significant and borderline status again
  # after quantile normalization
  for(type in c("totrna", "ribo", "Prot")){
    var1 = paste0("sig",type,"OnTP",i)
    var2 = paste0("borderlineZero",type,"OnTP",i)
    
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

# converting to foldchanges to rank
# but using the same colors
alldfsRank = list()
for(i in names(alldfs)){
  alldfsRank[[i]] = alldfs[[i]] %>% 
    mutate(mRNA = rank(mRNA),
           protein = rank(protein),
           RPF = rank(RPF))
}

# computing and plotting correlation for
# highest and lowest ranks
corPlot = function(df, type){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x) ; quad = paste0(type,"-quad")}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x); quad = paste0(type,"-quad")}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x); quad = paste0(type,"-quad")}
  
  cortop = tibble()
  
  dfhead = df %>%
    filter(get(quad) == "Q1") %>% 
    arrange(desc(get(y)))
  
  for(i in 1:dim(dfhead)[1]){
    cortop[i,"idx"] = i
    cortop[i,"value"] = cor(y=dfhead[1:i,y], x=dfhead[1:i,x], method = "pearson") %>% as.numeric()
    cortop[i,"pos"] = "up"
  }
  
  corbottom = tibble()
  
  dfhead = df %>%
    filter(get(quad) == "Q3") %>% 
    arrange(get(y))
  
  for(i in 1:dim(dfhead)[1]){
    corbottom[i,"idx"] = i
    corbottom[i,"value"] = cor(y=dfhead[1:i,y], x=dfhead[1:i,x], method = "pearson") %>% as.numeric()
    corbottom[i,"pos"] = "down"
  }
  
  cordf = bind_rows(corbottom,cortop) %>% drop_na()
  
  plot = ggplot(cordf, aes(x=idx,y=value,colour=pos)) +
                  geom_line(show.legend = F) +
    scale_color_manual(values = c("up"="#E15759", "down"="#4E79A7")) +
    ylim(c(-1,1))
  
  return(plot)
}

# plotting rank correlation trends
types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")
rankCorPlots = list()
for(i in names(alldfsRank)){
  for(k in types){
    rankCorPlots[[i]][[k]] = corPlot(alldfs[[i]], k)}
}
ggarrange(plotlist=c(rankCorPlots$TP4), nrow = 1, ncol = 3)

# including only those genes that have significant protein
# change in at least one timepoint
alldfsProtSigAtLeastOneTP = list()
for(i in names(alldfs)){
  filt = alldfs[[i]][,"locus_tag"] %>% unlist() %>% unname() %in% sigAtLeastOneTP
  alldfsProtSigAtLeastOneTP[[i]] = alldfs[[i]][filt,]
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

# preparing alldfs object
fulldf = list()
for(i in names(alldfs)){
  fulldf[[i]] = alldfs[[i]] %>% rename_at(vars(starts_with("sigProt")),
                                          list(~sub("TP[1-4]$", "", .)))
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

# plotting general tendencies
plotGeneral = function(df, type, xLim, yLim){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x)}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x)}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x)}
  
  model = lm(formula = form, data = df)
  cortext = paste0("R^2 = ", summary(model)$r.squared %>% round(digits = 3))
  p = pf(summary(model)$fstatistic[1],
         summary(model)$fstatistic[2],
         summary(model)$fstatistic[3],
         lower.tail = FALSE) %>% 
    formatC(format = "e", digits = 2)
  pvaltext = paste0("p = ", p)
  
  grob = grobTree(textGrob(paste0(cortext, "; ", pvaltext),
                           gp=gpar(fontsize=8),
                           x=0.075, y=0.95, hjust=0))
  
  alp = 0.05
  sz = 2
  
  plot = ggplot(data = df) +
    geom_point2(aes(x=get(x), y=get(y)), alpha = alp, size = sz) +
    xlim(xLim) + ylim(yLim) + xlab(x) + ylab(y) +
#    geom_vline(xintercept = 0, size = 0.3) +
#    geom_hline(yintercept = 0, size = 0.3) +
    geom_point2(inherit.aes = F,
                data = df,
                aes(x=get(x), y=get(y)),
                alpha = alp, size = sz,
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

# # plotting general tendencies animation
# fulldf = bind_rows(fulldf)
# 
# # computing lims
# xLim = findXLims(fulldf)
# yLim = findYLims(fulldf)
# 
# # animate general tendencies
# plotGeneralAnimate = function(fulldf, type, xLim, yLim){
#   if(type == "RPF-mRNA"){
#     y = "RPF"; x = "mRNA"}
#   if(type == "protein-mRNA"){
#     y = "protein"; x = "mRNA"}
#   if(type == "protein-RPF"){
#     y = "protein"; x = "RPF"}
#   
#   plot = ggplot(data = fulldf, aes_string(x=x, y=y, color="sigProtOn", group="locus_tag")) +
#     geom_point(alpha = 0.25, show.legend = F) +
#     geom_smooth(inherit.aes = F,
#                 data = fulldf,
#                 mapping=aes_string(x=x, y=y),
#                 formula = y ~ x,
#                 method = "lm",
#                 color = "black") +
#     scale_colour_manual(values = c("yes"="red", "no"="grey")) +
#     xlim(xLim) + ylim(yLim)
#     labs(title = "Time Point {closest_state} vs. Time Point 1") +
#     transition_states(timePoint, wrap = T)
#   return(plot)
# }
# 
# # types of plot we should get
# types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")
# 
# # plots
# generalIntPlots = list()
# for(i in types){
#   generalIntPlots[[i]] = plotGeneralAnimate(fulldf, i, xLim, yLim)
#   anim_save(animation = generalIntPlots[[i]],
#             filename = paste0("generaIntPlot_",i,".gif"))
# }

# plot functions with sig dataset
plotLM = function(df, type, xLim, yLim, rankTransformation){
  if(type == "RPF-mRNA"){
    y = "RPF"; x = "mRNA"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}
  if(type == "protein-mRNA"){
    y = "protein"; x = "mRNA"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}
  if(type == "protein-RPF"){
    y = "protein"; x = "RPF"; form = paste0(y,"~",x); color = paste0(type, "-quad"); alpha = paste0(y,"_",x,"_","sig")}

  if(y == "RPF"){var1 = "sigriboOn"}
  if(y == "protein"){var1 = "sigProtOn"}
  if(x == "mRNA"){var2 = "sigtotrnaOn"}
  if(x == "RPF"){var2 = "sigriboOn"}
  
  varName1 = df %>%
    select(starts_with(var1)) %>% names() %>% as.symbol()
  varName2 = df %>%
    select(starts_with(var2)) %>% names() %>% as.symbol()
  colsymb = color %>% as.symbol()
  dfSig = df %>% filter((!!varName1 == "yes" & !!varName2 == "yes") & (!!colsymb == "Q1" | !!colsymb == "Q3"))
#  dfNonSig =  df %>% filter(!!varName2 == "no")
  
#  alp = 0.25
  sz = 1.5
  
  var = color %>% as.symbol()
  model = list() ; cortext = list()
  p = list() ; pvaltext = list()
  grob = list() ; grobY = 0.92
  
  quadsIt = dfSig %>% select(!!var) %>% distinct() %>% unlist() %>% unname()
  quadsIt = quadsIt[quadsIt %in% paste0("Q", c(1,3))]
  
  for(i in quadsIt){
    model[[i]] = lm(formula = form, data = dfSig %>% filter(!!var == i))
    cortext[[i]] = paste0("R^2 = ", summary(model[[i]])$r.squared %>% round(digits = 3))
    p[[i]] = pf(summary(model[[i]])$fstatistic[1],
                summary(model[[i]])$fstatistic[2],
                summary(model[[i]])$fstatistic[3],
                lower.tail = FALSE) %>% 
      formatC(format = "e", digits = 2)
    pvaltext[[i]] = paste0("p = ", p[[i]])
    if(i == "Q1"){col = tab10col["Q1"]}
#    if(i == "Q2"){col = tab10col["Q2"]}
    if(i == "Q3"){col = tab10col["Q3"]}
#    if(i == "Q4"){col = tab10col["Q4"]}
    grob[[i]] = grobTree(textGrob(paste0(cortext[[i]], "; ", pvaltext[[i]]),
                                  gp=gpar(fontsize=8, col=col),
                                  x=0.075, y=grobY, hjust=0))
    grobY = grobY - 0.05
  }
  
  plot = ggplot(data = df, aes(x=get(x), y=get(y), colour=get(color))) +
    geom_point2(aes(alpha=get(alpha)), size = sz, show.legend = F) +
    xlim(xLim) + ylim(yLim) + xlab(x) + ylab(y) +
    # geom_vline(xintercept = 0, size = 0.3) +
    # geom_hline(yintercept = 0, size = 0.3) +
    # geom_point2(inherit.aes = F,
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
  
  varx = paste0(x,"Rank")
  df = df %>% arrange(get(x))
  dims = df %>% dim()
  neg = df[,x] < 0 ; neg = neg %>% sum() ; neg = (-neg+1):-1
  pos = 1:((neg %>% length() + 1):dims[1] %>% length())
  vec = c(neg,pos)
  df = df %>% mutate(!!varx := vec)
  minvarx = vec %>% min(); maxvarx = vec %>% max()

  vary = paste0(y,"Rank")
  df = df %>% arrange(get(y))
  dims = df %>% dim()
  neg = df[,y] < 0 ; neg = neg %>% sum() ; neg = (-neg+1):-1
  pos = 1:((neg %>% length() + 1):dims[1] %>% length())
  vec = c(neg,pos)
  df = df %>% mutate(!!vary := vec)
  minvary = vec %>% min(); maxvary = vec %>% max()
  
  rankplot = ggplot(data = df, aes(x=get(varx), y=get(vary), colour=get(color))) +
    geom_point2(aes(alpha=get(alpha)), size = sz, show.legend = F) +
    xlab(x) + ylab(y) + xlim(minvarx,maxvarx) + ylim(minvary,maxvary) +
    scale_colour_manual(values = tab10col) +
    scale_fill_manual(values = tab10col) +
    scale_alpha_manual(values = c("bold"=0.8, "mid"= 0.3, "faint"=0.1))
  
  if(rankTransformation == "yes"){
    return(rankplot)
  }else{
    return(plot)
  }
}

# types of plot we should get
types = c("RPF-mRNA", "protein-mRNA", "protein-RPF")

######## linear models and plots for full dataset
lmplots = list()
for(i in names(alldfs)){
  xLim = findXLims(alldfs[[i]])
  yLim = findYLims(alldfs[[i]])
  for(k in types){
    lmplots[[i]][[k]] = plotLM(alldfs[[i]], k, xLim, yLim, "no")
  }
}

# arranging plots (per time point)
p = list()
for(i in names(alldfsProtSigAtLeastOneTP)){
  p[[i]] = ggarrange(plotlist = c(lmplots[[i]]),
                nrow = 1, ncol = 3)
}

######## rankplots for full dataset
rankPlots = list()
for(i in names(alldfs)){
  xLim = findXLims(alldfs[[i]])
  yLim = findYLims(alldfs[[i]])
  for(k in types){
    rankPlots[[i]][[k]] = plotLM(alldfs[[i]], k, xLim, yLim, "yes")}
}
ggarrange(plotlist = c(rankPlots$TP4), nrow = 1, ncol = 3)

# arranging plots (per time point)
p = list()
for(i in names(alldfsProtSigAtLeastOneTP)){
  p[[i]] = ggarrange(plotlist = c(rankPlots[[i]]),
                     nrow = 1, ncol = 3)
}

# arranging plots in a panel
# ggarrange(plotlist = list(pmrnarpf$TP2, pmrna$TP2,prpf$TP2), ncol = 3, nrow = 1, labels = "TP2 vs TP1", legend = "none")
# ggarrange(plotlist = list(pmrnarpf$TP3,pmrna$TP3,prpf$TP3), ncol = 3, nrow = 1, labels = "TP3 vs TP1", legend = "none")

# arranging plots (full panel)
ggarrange(plotlist = c(lmplots$TP2,
                       lmplots$TP3,
                       lmplots$TP4),
          nrow = 3, ncol = 3)

# saving last object
#save(alldfsProtSigAtLeastOneTP, file = "alldfsProtSigAtLeastOneTP.RData")
