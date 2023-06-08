# PC3 cr 11 trees BIC
library(fs)
library(l1ou)
library(purrr)
library(geiger)
library(dplyr)
library(stringr)
source("2021_conv/R/modifiedFnsSummarize.R")

# all RDS files
resRDSs <- dir_ls("2021_conv/out/",regexp = "_bic")

# PC3
PC3res <- resRDSs[str_detect(resRDSs,"/emodPC3_cr")]
PC3res
PC3_l1oumods <-map(PC3res,readRDS)
PC3_l1oumods
fp0 <- PC3_l1oumods[1:11]
fp1 <- flatten(PC3_l1oumods[12])
fp2 <- flatten(PC3_l1oumods[13])

allpc3 <- c(fp0,fp1,fp2)

PC3_conv_regs <- map(allpc3,estimate_convergent_regimes)
saveRDS(PC3_conv_regs,"2021_conv/out/PC3_conv_regs_cr-all100trees.rds")

# condense convergent clades
convtipsPC3 <- map(PC3_conv_regs,summarize.lasso.convergence.output)
convtipslongPC3 <- convtipsPC3 %>% map_df(bind_rows,.id = "tree")

justcon_PC3 <- convtipsPC3 %>% map(~filter(.x,conv_reg!="noshift"))
justcon_PC3long <- justcon_PC3 %>% map_df (bind_rows,.id="tree")


# all species
allsps_conv <- convtipslongPC3 %>% distinct(sp) %>% pull(sp)

# 'basal' species across all trees
treesgp <- convtipslongPC3 %>% group_by(tree) %>% n_groups()
noShift_taxa <- convtipslongPC3 %>% group_by(tree) %>% filter(conv_reg=="noshift") %>% 
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
  convtipsPC3 %>% 
    map(~get_conv_taxa(.x,focal_sp)) %>% 
    map_df('convergent_sps') %>% add_count(sp)
}

PC3allbros <- map(convtaxa,count_bros) %>% map_df(bind_rows) 
PC3allbros <- PC3allbros %>% distinct()
PC3allbros
PC3allbros %>% as_tbl_graph()

PC3allbros %>% rename(to=1,from=2,weight=3) %>% 
  as_tbl_graph()-> adj_grph


PC3allbros %>% 
  filter(n>4) %>% 
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
