# Process 1lou models into sets of edges for graph visualization
library(fs)
library(l1ou)
library(purrr)
library(geiger)
library(dplyr)
library(stringr)
source("2021_conv/R/modifiedFnsSummarize.R")
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(ggrepel)
library(ggtree)

# Load PCs cr convergent regimes----
# PC1cr 
PC1_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC1_convRegs_crAll.rds")
PC2_conv_regs_cr  <- readRDS("2021_conv/out/processed_l1ouModels/PC2_convRegs_crAll.rds")
PC3_conv_regs_cr  <- readRDS("2021_conv/out/processed_l1ouModels/PC3_convRegs_crAll.rds")
# Multivariate cr
conv_multcr <- readRDS("2021_conv/out/processed_l1ouModels/mult_convRegs_crAll.rds")
# ext convergent regimes
PC1_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC1_convRegs_extAll.rds")
PC2_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC2_convRegs_extAll.rds")
PC3_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC3_convRegs_extAll.rds")
# multivariate ext conv regs
conv_multext <- readRDS("2021_conv/out/processed_l1ouModels/mult_convRegs_extAll.rds")

# all measurement PCs
conv_multALL <- readRDS("2021_conv/out/processed_l1ouModels/ALLPCs_convRegs_All.rds")

# using custom function to collapse convergent clades
PC1allcr <- summarise_convregs(PC1_conv_regs_cr)
PC1allext <- summarise_convregs(PC1_conv_regs_ext)
PC2allcr <- summarise_convregs(PC2_conv_regs_cr)
PC2allext <- summarise_convregs(PC2_conv_regs_ext)
PC3allcr <- summarise_convregs(PC3_conv_regs_cr)
PC3allext <- summarise_convregs(PC3_conv_regs_ext)
multAllcr <- summarise_convregs(conv_multcr)
multAllext <- summarise_convregs(conv_multext)
allPCsMult <- summarise_convregs(conv_multALL)


write.csv(PC1allcr,"2021_conv/out/processed_edges/crPC1resEdges.csv")
write.csv(PC2allcr,"2021_conv/out/processed_edges/crPC2resEdges.csv")
write.csv(PC3allcr,"2021_conv/out/processed_edges/crPC3resEdges.csv")
write.csv(multAllcr,"2021_conv/out/processed_edges/crMultresEdges.csv")
write.csv(PC1allext,"2021_conv/out/processed_edges/extPC1resEdges.csv")
write.csv(PC2allext,"2021_conv/out/processed_edges/extPC2resEdges.csv")
write.csv(PC3allext,"2021_conv/out/processed_edges/extPC3resEdges.csv")
write.csv(multAllext,"2021_conv/out/processed_edges/extMultresEdges.csv")
write.csv(allPCsMult,"2021_conv/out/processed_edges/allPCsMultresEdges.csv")

# redo
PC1allcr <- read.csv("2021_conv/out/processed_edges/crPC1resEdges.csv") %>% select(-1)
PC2allcr <- read.csv("2021_conv/out/processed_edges/crPC2resEdges.csv") %>% select(-1)
PC3allcr <- read.csv("2021_conv/out/processed_edges/crPC3resEdges.csv") %>% select(-1)
multAllcr <- read.csv("2021_conv/out/processed_edges/crMultresEdges.csv") %>% select(-1)
PC1allext <- read.csv("2021_conv/out/processed_edges/extPC1resEdges.csv") %>% select(-1)
PC2allext <- read.csv("2021_conv/out/processed_edges/extPC2resEdges.csv") %>% select(-1)
PC3allext <- read.csv("2021_conv/out/processed_edges/extPC3resEdges.csv") %>% select(-1)
multAllext <- read.csv("2021_conv/out/processed_edges/extMultresEdges.csv") %>% select(-1)
allPCsMult <- read.csv("2021_conv/out/processed_edges/allPCsMultresEdges.csv") %>% select(-1)

ggplot(PC1allcr)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(PC1allext)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(PC2allcr)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(PC2allext)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(PC3allcr)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(PC1allext)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(multAllcr)+
  geom_histogram(aes(x=n),binwidth = 1)
ggplot(multAllext)+
  geom_histogram(aes(x=n),binwidth = 1)

# filter to get concensus (CR) ----
PC1allcr %>%
  filter(n > quantile(PC1allcr$n,0.90)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphPC1
adj_grphPC1

PC1allcr_regimeMembership <- clusters(adj_grphPC1)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
components(adj_grphPC1)

grpcomps <- split(V(adj_grphPC1), components(adj_grphPC1)$membership)
gp1 <- subgraph(adj_grphPC1,grpcomps[[1]]) |> as_tbl_graph() 
plot(gp1)
cliques(gp1)

plot(gp1)
edgesgp1 <- gp1 %>% get.edgelist() |> as_tibble()

gp2 <- subgraph(adj_grphPC1,grpcomps[[2]]) |> as_tbl_graph() 

PC2allcr %>%
  filter(n > quantile(PC2allcr$n,0.9)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphPC2
PC2allcr_regimeMembership <- clusters(adj_grphPC2)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2)
PC2allcr_regimeMembership

PC3allcr %>%
  filter(n > quantile(PC3allcr$n,0.9)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphPC3
PC3allcr_regimeMembership <- clusters(adj_grphPC3)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2)
PC3allcr_regimeMembership

multAllcr %>%
  filter(n > quantile(multAllcr$n,0.90)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphMult

multAllcr_regimeMembership <- clusters(adj_grphMult)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2)
multAllcr_regimeMembership

# filter for concensus (ext) ----

adj_grphPC1ext
PC1allext_regimeMembership <- clusters(adj_grphPC1ext)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
components(adj_grphPC1ext)

PC2allext %>%
  filter(n >= quantile(PC2allext$n,0.9)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphPC2ext
adj_grphPC2ext
PC2allext_regimeMembership <- clusters(adj_grphPC2ext)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
components(adj_grphPC2ext)



PC2allcr %>%
  filter(n > quantile(PC2allcr$n,0.9)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphPC2
PC2allcr_regimeMembership <- clusters(adj_grphPC2)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2)
PC2allcr_regimeMembership

PC3allext %>%
  filter(n >= quantile(PC3allext$n,0.9)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphPC3ext
PC3allext_regimeMembership <- clusters(adj_grphPC3ext)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2)
PC3allext_regimeMembership

multAllext %>%
  filter(n > quantile(multAllext$n,0.90)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphMultext
multAllext_regimeMembership <- clusters(adj_grphMultext)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2)
multAllext_regimeMembership

# whole body
allPCsMult %>%
  filter(n >= quantile(allPCsMult$n,0.90)) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph() -> adj_grphMultAll



ggraph(adj_grphMultAll, layout = "fr") +
  geom_edge_link(aes(color = factor(treeID))) +
  geom_node_point(size = 2) +
  theme_graph() +
  geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 3)


# adj_grphMultAll %>%
#     activate(edges) |>
#   mutate(group = tidygraph::group_biconnected_component()) |> as_tibble()
# adj_grphMultAll %>%
# #  activate(edges) |>
#   mutate(group = tidygraph::group_infomap()) |> 
#  # activate(nodes) |> 
# ggraph(layout = "stress", weights = weight) +
#   geom_edge_fan(aes(color = factor(weight), alpha = factor(weight))) +
#   geom_node_point(size = 2,aes(color=factor(group))) +
#   theme_graph() +
#   geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 3)

components(adj_grphMultAll)

# PC intersections----  
# CR
PC1pairs <- PC1allcr %>%
  filter(n > quantile(PC1allcr$n,0.9)) |> select(-n)
PC1pairs <- PC1pairs[!duplicated(t(apply(PC1pairs, 1, sort))),]
PC1pairs

PC2pairs <- PC2allcr %>%
  filter(n > quantile(PC2allcr$n,0.9)) |> select(-n)
PC2pairs <- PC2pairs[!duplicated(t(apply(PC2pairs, 1, sort))),]

PC3pairs <- PC3allcr %>%
  filter(n > quantile(PC3allcr$n,0.9)) |> select(-n)
PC3pairs <- PC3pairs[!duplicated(t(apply(PC3pairs, 1, sort))),]

Multpairs <- multAllcr %>%
  filter(n > quantile(multAllcr$n,0.9)) |> select(-n)
Multpairs<- Multpairs[!duplicated(t(apply(Multpairs, 1, sort))),]
Multpairs

# CR
PC12int <- intersect(PC1pairs,PC2pairs)
PC12int
PC12int[!duplicated(t(apply(PC12int, 1, sort))),]
PC12int |> as_tbl_graph() -> adj_grphPC12int
PC13int <- intersect(PC1pairs,PC3pairs)
PC13int
PC23int <- intersect(PC2pairs,PC3pairs)
PC23int
intersect(PC1pairs,Multpairs)
intersect(PC2pairs,Multpairs)
intersect(PC3pairs,Multpairs)

ggraph(adj_grphPC12int, layout = "stress") +
  geom_edge_link() +
  geom_node_point(size = 1) +
  theme_graph() +
  geom_text_repel(aes(x,y,label = name),fontface="italic",color='black',bg.color='white',bg.r=0.15)

# EXT
PC1pairsext <- PC1allext %>%
  filter(n > quantile(PC1allext$n,0.9)) |> select(-n)
PC1pairsext <- PC1pairsext[!duplicated(t(apply(PC1pairsext, 1, sort))),]
PC1pairsext

PC2pairsext <- PC2allext %>%
  filter(n >= quantile(PC2allext$n,0.9)) |> select(-n)
PC2pairsext <- PC2pairsext[!duplicated(t(apply(PC2pairsext, 1, sort))),]
PC2pairsext

PC3pairsext <- PC3allext %>%
  filter(n >= quantile(PC3allext$n,0.9)) |> select(-n)
PC3pairsext <- PC3pairsext[!duplicated(t(apply(PC3pairsext, 1, sort))),]
PC3pairsext

Multpairsext <- multAllext |> 
  filter(n >= quantile(multAllext$n,0.9)) |> select(-n)
Multpairsext <- Multpairsext[!duplicated(t(apply(Multpairsext, 1, sort))),]
Multpairsext

PC3pairsext
PC12intext <- intersect(PC1pairsext,PC2pairsext)
PC12intext
PC13intext <- intersect(PC1pairsext,PC3pairsext)
PC13intext
PC23intext <- intersect(PC2pairsext,PC3pairsext)
PC23intext
intersect(PC1pairsext,Multpairsext)
intersect(PC2pairsext,Multpairsext)
intersect(PC3pairsext,Multpairsext)


# CR vs EXT
intersect(PC1pairs,PC1pairsext)
intersect(PC1pairs,PC2pairsext)
intersect(PC1pairs,PC3pairsext)
intersect(PC1pairs,Multpairsext)
intersect(PC2pairs,PC1pairsext)
intersect(PC2pairs,PC2pairsext)
intersect(PC2pairs,PC3pairsext)
intersect(PC2pairs,Multpairsext)
intersect(PC3pairs,PC1pairsext)
intersect(PC3pairs,PC2pairsext)
intersect(PC3pairs,PC3pairsext)
intersect(PC3pairs,Multpairsext)
####################################################################
# phylo
rodtrees <- read.nexus("2021_conv/phylo/rodent_trees265.nex")

# random tree
rodphly <- rodtrees[[34]]

# the ggtree groupOTU function seems to take lists, easy to do with split
groupInfoPC1cr <- 
  split(PC1allcr_regimeMembership$sp,PC1allcr_regimeMembership$regime)
# join Regs with the tree
PC1crPhylRegs <- groupOTU(rodphly, groupInfoPC1cr, overlap="abandon")
PC1crPhylRegs
codist <- cophenetic.phylo(PC1crPhylRegs)

xy <- t(combn(colnames(codist), 2))
data.frame(xy, dist=codist[xy])

edgesgp1 |> mutate(across(everything(),~str_replace(.," ","_"))) |> 

ggtree(PC1crPhylRegs,layout = "fan",aes(color=group))+
  scale_color_manual(values=c("black","#e7298a","#1f78b4"),
                     labels=c("1","2"),
                     breaks=c("1","2"))#+
  theme(legend.position="right")+xlim(0,1)

# MRCAs
# vector of tip lables 
tipstr <- PC1crPhylRegs$tip.label

# define function for getting MRCA node
getLineage_node <- function(regime){
    nodenum <- MRCA(PrimatesDat, tip=c(str_subset(tipstr,genusName)))
    outt <- tibble(genus=genusName,node=nodenum)
    return(outt)
}
  
  # table with genera and node numbers
  genNodes <- map(prGenera,getgenusNode) %>% bind_rows
  
  # loop to draw the labels
  for(j in 1:nrow(genNodes)){
    #Then add each clade label
    gheat <- gheat + geom_cladelabel(node=genNodes$node[j], label=genNodes$genus[j], 
                                     barsize = 1, offset.text = 0.02, offset =0.06,angle = "auto",
                                     family = "Roboto", fontsize=2.4)
  }
  
  gheat
  
  
  
  
  #### 
  grpcomps <- split(V(adj_grphPC1), components(adj_grphPC1)$membership)
  gp1 <- subgraph(adj_grphPC1,grpcomps[[1]]) |> as_tbl_graph() 
  plot(gp1)
  cliques(gp1)
  
  plot(gp1)
  edgesgp1 <- gp1 %>% get.edgelist() |> as_tibble()
  
  gp2 <- subgraph(adj_grphPC1,grpcomps[[2]]) |> as_tbl_graph() 
  
  gp1 |> activate("edges") |> get.edgelist()|> as_tibble()
  gp2 |> activate("edges") |> get.edgelist()|> as_tibble()
  
  
  class(adj_grphPC1)
  grpcomps$`1`
  class(grpcomps[[1]])
  cliques(grpcomps[[1]])