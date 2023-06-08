# for gheatmap
library(purrr)
library(tidyr)
#library(ggh4x)
library(forcats)
library(aplot)

# special case for wb models
wbregsletters <- Multallwbtagged %>%
  pivot_longer(-treeID) %>% 
  rename(sp=value) %>% select(-name) %>% 
  distinct(sp,treeID) %>% 
  select(sp,treeID) 

# sp with two or more regs
wbdf <- 
wbregsletters %>% mutate(PC="MultALLPCs",
                    regime=as.numeric(chartr("abcdef","123456",treeID))) %>% 
  select(-treeID) 

regslist <- list(
PC1regscr, 
PC2regscr, 
PC3regscr, 
Multregscr, 
PC1regsext,
PC2regsext, 
PC3regsext, 
Multregsext) 


allregsdf <- purrr::map_df(regslist,1)
allregsdf
allregsdfwide <- allregsdf %>% select(-mtype) %>% 
  pivot_wider(names_from = PC,values_from = regime)

# Phylogeny figures -------------------------------------------------------
# phylo

allregmat <- allregsdfwide %>% mutate(across(where(is.numeric),as.factor)) 
allregmat
allregmat <- allregmat %>% as.data.frame()
row.names(allregmat) <- allregmat$sp
allregmat$sp <- NULL
#rectphyl <- ggtree(rodphyl,layout = "rectangular")
#rectaggr <- pcfig %<+% allregsdfwide
#rectaggr
#gheatmap(rectaggr,allregmat, width = 1,offset=8)

#allreg4tile <- allregmat %>% tibble::rowid_to_column()

allregsdf <- 
allregsdf %>% 
  mutate(grplab=if_else(mtype=="cranial","craniodental & mandibular",mtype))

# cr and external, get all tips
treetaxadf <- tibble(treetips=get_taxa_name(pcfig))
allregsdf <- left_join(treetaxadf,allregsdf,by=c("treetips"="sp"))
allregsdf$sp <- factor(allregsdf$treetips,levels = get_taxa_name(pcfig))
allregsFilt <- 
  allregsdf %>% filter(!is.na(mtype))
allregsFilt$PC <- factor(allregsFilt$PC,levels = c("PC1cr","PC2cr","PC3cr","Multcr",
                                 "PC1ext","PC2ext","PC3ext","Multext"))

# prepare wb data for its own tile chart
allwbregs <- left_join(treetaxadf,wbdf,by=c("treetips"="sp"))
allwbregs$sp <- factor(allwbregs$treetips,levels = get_taxa_name(pcfig))
allwbregs$PC <- factor(allwbregs$PC,levels = c("MultALLPCs"))
allwbregs$regime <- factor(allwbregs$regime,levels = c("1","2","3","4","5","6"))
allwbregsFilt <- allwbregs %>% filter(!is.na(regime))
allwbregsFilt$facetlab <- "whole body"

# plot cr and ext 
tiles <- 
ggplot(allregsFilt)+
  geom_tile(aes(y=fct_rev(sp),x=PC,fill=factor(regime)),color="gray")+
  facet_wrap(~grplab,
             scales = "free",drop=TRUE,
             labeller = labeller(grplab = label_wrap_gen(15)))+
  theme_minimal()+
     theme(
       strip.text = element_text(size=9,vjust = 0),
       text = element_text(family="Laksaman"),
         panel.spacing.x = unit(0.01, "lines"),
         axis.title.x = element_blank(), 
         axis.text.x = element_text(angle = 90),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         plot.margin = margin(0,0,0,-10),
         panel.grid.major.x = element_blank())+
   scale_y_discrete(drop=FALSE)+
   scale_x_discrete(labels=c("PC1","PC2","PC2","PCs 1-3"))+
   # coord_fixed(0.5)+
        scale_fill_scico_d(palette = "imola",guide="none")
tiles

# plot wb
wbtiles <- 
allwbregsFilt %>%
  ggplot(aes(x = PC, y = fct_rev(sp), fill = regime,group=regime)) +
  geom_tile(position = "dodge")+
  scale_y_discrete(drop=FALSE)+
  facet_wrap(~facetlab)+
  theme_minimal()+
  theme(
    strip.text = element_text(size=9),
    text = element_text(family="Laksaman"),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0,0,0,-10),
    panel.grid.major.x = element_blank())+
  scale_x_discrete(labels=c("PCs1-5"))+
  # coord_fixed(0.5)+
  scale_fill_scico_d(palette = "buda",guide="none")
wbtiles

# compose
tiles %>% insert_left(pcfig,width = 2) %>% insert_right(wbtiles,width = 0.5)

ggsave("2021_conv/figs/new/phyloDist.png",
       device = agg_png,
       height = 8,width = 6.7, units = "in",dpi = 400)
