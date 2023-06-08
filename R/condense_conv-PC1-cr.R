# PC1 all 100 trees BIC
library(fs)
#remotes::install_github("khabbazian/l1ou")
library(l1ou)
library(purrr)
library(geiger)
library(dplyr)
library(stringr)
source("2021_conv/R/modifiedFnsSummarize.R")
library(tidygraph)
library(ggraph)
library(graphlayouts)

# Pending conv. regimes ----
# # all RDS files
# resRDSs <- dir_ls("2021_conv/out/",regexp = "_bic")
# # PC1
# PC1res <- resRDSs[str_detect(resRDSs,"/emodPC1_cr")]
# PC1res <- PC1res[str_detect(PC1res,"trees",negate = TRUE)]
# PC1_l1oumods <- map(PC1res,readRDS)
# PC1_conv_regs <- map(PC1_l1oumods,estimate_convergent_regimes,nCores=4)
# saveRDS(object = PC1_conv_regs,file = "2021_conv/out/PC1_conv_regs_cr-11trees.rds")

# Load PC1 cr convergent regimes----
# releant RDS files
PC1convRDSs <- dir_ls("2021_conv/out/processed_l1ouModels/",regexp = "PC1_")
PC1convRDSs
PC1_conv_regs_cr <- readRDS(PC1convRDSs)


# condense convergent clades
convtipsPC1 <- map(PC1_conv_regs_cr,summarize.lasso.convergence.output)
convtipslongPC1 <- convtipsPC1 %>% map_df(bind_rows,.id = "tree")
justcon_PC1 <- convtipsPC1 %>% map(~filter(.x, conv_reg!="noshift"&regime_type=="conv_reg"))
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



# count members
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

hist(PC1allbros$n)

PC1allbros %>% 
  filter(n>7) %>% 
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
