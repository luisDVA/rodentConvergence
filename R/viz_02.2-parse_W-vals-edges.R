# Parse edges and W values for all the regimes
# prerequisite for adjacency plots
library(fs)
library(igraph)
library(purrr)
library(geiger)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(ggrepel)
library(ggtree)
library(scico)
library(patchwork)
library(ragg)
library(extrafont)
source("2021_conv/R/viz_02-helper_fns.R")

# Import processed edges data
PC1allcr <- read.csv("2021_conv/out/processed_edges/crPC1resEdges.csv") %>% select(-1)
PC2allcr <- read.csv("2021_conv/out/processed_edges/crPC2resEdges.csv") %>% select(-1)
PC3allcr <- read.csv("2021_conv/out/processed_edges/crPC3resEdges.csv") %>% select(-1)
Multallcr <- read.csv("2021_conv/out/processed_edges/crMultresEdges.csv") %>% select(-1)
PC1allext <- read.csv("2021_conv/out/processed_edges/extPC1resEdges.csv") %>% select(-1)
PC2allext <- read.csv("2021_conv/out/processed_edges/extPC2resEdges.csv") %>% select(-1)
PC3allext <- read.csv("2021_conv/out/processed_edges/extPC3resEdges.csv") %>% select(-1)
Multallext <- read.csv("2021_conv/out/processed_edges/extMultresEdges.csv") %>% select(-1)
Multallwbtagged <- read.csv("2021_conv/out/processed_edges/allPCsMultresEdgesGrouped.csv") 

# filter to get concensus (CR) ----
PC1regscr <- to_regComponents(PC1allcr,">","PC1cr","cranial")
PC2regscr <- to_regComponents(PC2allcr,">","PC2cr","cranial")
PC3regscr <- to_regComponents(PC3allcr,">","PC3cr","cranial")
Multregscr <- to_regComponents(Multallcr,">","Multcr","cranial")
PC1regsext <- to_regComponents(PC1allext,">","PC1ext","external")
PC2regsext <- to_regComponents(PC2allext,">=","PC2ext","external")
PC3regsext <- to_regComponents(PC3allext,">=","PC3ext","external")
Multregsext <- to_regComponents(Multallext,">","Multext","external")
Multregswb <- to_regComponentswb(Multallwbtagged,"MultALLPCs","wb")

# read and summarize w values
wresults <- dir_ls("2021_conv/out/processed_regWs/",glob = "*.csv")
allW_vals <- map_df(wresults,read_csv)
allW_vals <- allW_vals %>% mutate(mtype=
                                    case_when(str_detect(PC,"cr")~"cranial",
                                              str_detect(PC,"ext")~"external",
                                              TRUE~"wb"))
regWs_summary <- allW_vals %>% group_by(reg,PC) %>% 
  summarize(w=mean(w)) %>%
  ungroup()

# join regimes (group components)
PC1regscrGrph <- makeWgraph(PC1regscr,"PC1cr")
PC2regscrGrph <- makeWgraph(PC2regscr,"PC2cr")
PC3regscrGrph <- makeWgraph(PC3regscr,"PC3cr")
MultregscrGrph <- makeWgraph(Multregscr,"Multcr")
PC1regsextGrph <- makeWgraph(PC1regsext,"PC1ext")
PC2regsextGrph <- makeWgraph(PC2regsext,"PC2ext")
PC3regsextGrph <- makeWgraph(PC3regsext,"PC3ext")
MultregsextGrph <- makeWgraph(Multregsext,"Multext")
wbregsGrph <- makeWgraphwb(Multregswb,"MultALLPCs")


