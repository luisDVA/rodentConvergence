# Cranial, all 100 trees BIC, model data
library(fs)
library(l1ou)
library(purrr)
library(geiger)
library(dplyr)
library(stringr)
source("2021_conv/R/modifiedFnsSummarize.R")
library(ggplot2)
library(ggridges)
library(ggraph)
library(tidygraph)
library(tidyr)

# Load PC1 cr convergent regimes----
# releant RDS files
PC1_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC1_convRegs_crAll.rds")

# all shifts 
# Exploratory
PC1_conv_regs_cr %>% map_int('nShifts') %>% median
PC1_conv_regs_cr %>% map_int('nShifts') %>% sd %>% round(2)
PC1_conv_regs_cr %>% map_int('nShifts')  %>% min
PC1_conv_regs_cr %>% map_int('nShifts')  %>% max


PC1crRegimes <- PC1_conv_regs_cr %>% map_df(regimes_summary) %>% mutate(mod="PC1")

# Load PC2 cr convergent regimes----
PC2_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC2_convRegs_crAll.rds")
PC2crRegimes <- PC2_conv_regs_cr %>% map_df(regimes_summary) %>% mutate(mod="PC2")

# Load PC3 cr convergent regimes----
PC3_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC3_convRegs_crAll.rds")
PC3crRegimes <- PC3_conv_regs_cr %>% map_df(regimes_summary)%>% mutate(mod="PC3")

# multivariate cr conv regs
allPCs_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/mult_convRegs_crAll.rds")
allPCscrRegimes <- allPCs_conv_regs_cr %>% map_df(regimes_summary)%>% mutate(mod="multivariate")


regresultsCR <- bind_rows(PC1crRegimes,PC2crRegimes,PC3crRegimes,allPCscrRegimes)

regs_summary <- 
regresultsCR %>% group_by(mod) %>% summarise(nshifts=median(shifts),
                                             minshifts=min(shifts),
                                             maxshifts=max(shifts),
                                             nregimes=median(regimes),
                                             minregimes=min(regimes),
                                             maxregimes=max(regimes),
                                             nconv_regimes=median(conv_regimes),
                                             minconvregs=min(conv_regimes),
                                             maxconvregs=max(conv_regimes)) # %>% 
  

regs_summary %>% mutate(nshifts=paste0(nshifts," (",minshifts,"-",maxshifts,")"),
                        nregimes=paste0(nregimes," (",minregimes,"-",maxregimes,")"),
                        nconv_regimes=paste0(nconv_regimes," (",minconvregs,"-",maxconvregs,")")) %>% 
  select(mod,nshifts,nregimes,nconv_regimes) %>% gt::gt()

# all shifts 
PC1_conv_regs_cr %>% map_int('nShifts') %>% median
PC1_conv_regs_cr %>% map_int('nShifts') %>% sd %>% round(2)
PC1_conv_regs_cr %>% map_int('nShifts')  %>% min
PC1_conv_regs_cr %>% map_int('nShifts')  %>% max

regresultsCR

regresultsCR %>% pivot_longer(-mod) %>% 
  ggplot(aes(y=name,x=value,height = stat(density)))+
  geom_density_ridges(stat="binline",bins=8,scale=0.89)+
  facet_grid(~mod)





# condense convergent clades
convtipsPC1 <- map(PC1_conv_regs_cr,summarize.lasso.convergence.output)
convtipslongPC1 <- convtipsPC1 %>% map_df(bind_rows,.id = "tree")

justcon_PC1 <- convtipsPC1 %>% map(~filter(.x,conv_reg!="noshift"))
justcon_PC1long <- justcon_PC1 %>% map_df (bind_rows,.id="tree")


# all species
allsps_conv <- convtipslongPC1 %>% distinct(sp) %>% pull(sp)

# 'basal' species across all trees
treesgp <- convtipslongPC1 %>% group_by(tree) %>% n_groups()
noShift_taxa <- convtipslongPC1 %>% group_by(tree) %>% filter(conv_reg=="noshift") %>% 
  ungroup() %>% count(sp) %>% filter(n==treesgp) %>% pull(sp)

convtaxa <- setdiff(allsps_conv,noShift_taxa)


# get all convergent taxa for a focal sp, remove sps from same regime
get_conv_taxa <- function(modeldf,species){
  mbros <- modeldf %>% filter(stringr::str_detect(sp,species))
  n_basal <- nrow(mbros[mbros$conv_reg=="noshift",])
  mbros_convOnly <- mbros %>%  filter(conv_reg!="noshift") 
  conv_regimeN <- pull(mbros_convOnly,conv_reg)
  og_regimeN <- pull(mbros_convOnly,og_reg)
  convBros <-  modeldf %>% filter(conv_reg %in% conv_regimeN)
  conv_sps <- 
    convBros %>% filter(sp!=species) %>% 
    filter(!og_reg %in% og_regimeN) %>% 
    mutate(focal=species,.before=1) %>% 
    select(focal,sp) %>% distinct()
  return(list(noshift=n_basal,convergent_sps=conv_sps))
}



#
count_bros <- function(focal_sp){
  convtipsPC1 %>% 
    map(~get_conv_taxa(.x,focal_sp)) %>% 
    map_df('convergent_sps') %>% add_count(sp)
}

PC1allbros <- map(convtaxa,count_bros) %>% map_df(bind_rows) 
PC1allbros <- PC1allbros %>% distinct()
PC1allbros
PC1allbros %>% as_tbl_graph()

PC1allbros %>% rename(to=1,from=2,weight=3) %>% 
  as_tbl_graph()-> adj_grph

adj_grph
hist(PC1allbros$n)

PC1allbros %>% 
  filter(n>50) %>% 
  rename(to=1,from=2,weight=3) %>% 
  mutate(weight=weight)%>% 
  as_tbl_graph()-> adj_grph

adj_grph

ggraph(adj_grph,layout="kk",weights=weight)+
  geom_edge_link(aes(color=factor(weight),alpha=weight))+
  geom_node_point(size = 2)+
  theme_graph()+
  geom_node_text(aes(label = name), repel = TRUE,color="black",size=3)

ggraph(adj_grph, layout = "fr")+
  geom_edge_link(colour = "grey")+
  geom_node_point(size = 2)+
  theme_graph()

gtid %>% 
  activate(nodes) %>%
  mutate(Cluster = as.factor(group_components())) %>% 
  ggraph(layout = "kk") + 
  geom_edge_link(width = 0.1, colour = "lightgray") +
  geom_node_point(aes(colour = Cluster), size = 4) +
  scale_color_scico_d(palette = "davos")+
  theme_graph(background = "black",text_colour = "white")
