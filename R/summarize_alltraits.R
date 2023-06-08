# whole body model (special case)
# summarize convergent species pairs
library(fs)
library(l1ou)
library(purrr)
library(geiger)
library(dplyr)
library(stringr)
source("2021_conv/R/slat2019fns.R")
source("2021_conv/R/modifiedFnsSummarize.R")
library(ggplot2)
library(ggridges)
library(tidyr)
library(ggraph)
library(tidygraph)

# releant RDS files
# Load allPCS  convergent regimes--
ALLPCS_conv_regs <- readRDS("2021_conv/out/processed_l1ouModels/ALLPCs_convRegs_All.rds")
# summarize
ALLPCSRegimes <- ALLPCS_conv_regs %>% map_df(regimes_summary) %>% mutate(mod="ALLPCS")

# report shifts
regs_summary <- 
  ALLPCSRegimes %>%  summarise(nshifts=median(shifts),
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
  mutate(model="multivariate\n (all measurements)") %>% 
  select(model,nshifts,nregimes,nconv_regimes) %>% gt::gt()


# for individual convergence reporting
conv_vec <- ALLPCSRegimes |> tibble::rowid_to_column() |> filter(conv_regimes>0) |> pull(rowid)
convResAll <- ALLPCS_conv_regs[conv_vec]

# iterate each tree with conv. regs.
summarise_convregs(convResAll[1]) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphAllt1
summarise_convregs(convResAll[2]) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphAllt2
summarise_convregs(convResAll[3]) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphAllt3
summarise_convregs(convResAll[4]) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphAllt4
summarise_convregs(convResAll[5]) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphAllt5
summarise_convregs(convResAll[6]) %>%
  rename(to = 1, from = 2, weight = 3) %>%
  mutate(weight = weight) %>%
  mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
  as_tbl_graph(directed=FALSE) -> adj_grphAllt6

# intersect (for reporting)
T1pairs <- summarise_convregs(convResAll[1]) %>%  select(-n)
T1pairs <- T1pairs[!duplicated(t(apply(T1pairs, 1, sort))),]
T1pairs$treeID <- "a"
T2pairs <- summarise_convregs(convResAll[2]) %>%  select(-n)
T2pairs <- T1pairs[!duplicated(t(apply(T2pairs, 1, sort))),]
T2pairs$treeID <- "b"
T3pairs <- summarise_convregs(convResAll[3]) %>%  select(-n)
T3pairs <- T3pairs[!duplicated(t(apply(T3pairs, 1, sort))),]
T3pairs$treeID <- "c"
T4pairs <- summarise_convregs(convResAll[4]) %>%  select(-n)
T4pairs <- T4pairs[!duplicated(t(apply(T4pairs, 1, sort))),]
T4pairs$treeID <- "d"
T5pairs <- summarise_convregs(convResAll[5]) %>%  select(-n)
T5pairs <- T5pairs[!duplicated(t(apply(T5pairs, 1, sort))),]
T5pairs$treeID <- "e"
T6pairs <- summarise_convregs(convResAll[6]) %>%  select(-n)
T6pairs <- T6pairs[!duplicated(t(apply(T6pairs, 1, sort))),]
T6pairs$treeID <- "f"

# check for shared species pairs
intersect(T1pairs,T2pairs)
intersect(T1pairs,T3pairs)
intersect(T1pairs,T4pairs)
intersect(T1pairs,T5pairs)
intersect(T1pairs,T6pairs)
intersect(T2pairs,T3pairs)
intersect(T2pairs,T4pairs)
intersect(T2pairs,T5pairs)
intersect(T2pairs,T6pairs)
intersect(T3pairs,T4pairs)
intersect(T3pairs,T5pairs)
intersect(T3pairs,T6pairs)
intersect(T4pairs,T5pairs)
intersect(T5pairs,T6pairs)

# Get a list of all edges in all graphs
allTreesMult <- bind_rows(T1pairs,T2pairs,T3pairs,
                        T4pairs,T5pairs,T6pairs)


allTreesMult
readr::write_csv(allTreesMult,file = "2021_conv/out/processed_edges/allPCsMultresEdgesGrouped.csv")
allTreesMult |> as_tbl_graph() -> adj_grphMultAll
ggraph(adj_grphMultAll, layout = "fr") +
  geom_edge_link(aes(color = factor(treeID))) +
  geom_node_point(size = 2) +
  theme_graph() +
  geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 3)

intersect(igraph::as_data_frame(adj_grphAllt2),igraph::as_data_frame(adj_grphAllt5))

#
(union(adj_grphAllt1,adj_grphAllt2)) #,adj_grphAllt3,
#  adj_grphAllt4,adj_grphAllt5,adj_grphAllt6))

# Make a graph of all of the edges including overlap
uniongraph <- graph(do.call(rbind, alledges))
uniongraph

 resultgraph <-
  graph.adjacency(get.adjacency(uniongraph) > 0.1 * 6)
resultgraph
ggraph(uniongraph, layout = "nicely") +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph() +
  geom_node_text(aes(label = name), repel = TRUE, color = "black", size = 3)

intersect(igraph::as_data_frame(adj_grphPC1),igraph::as_data_frame(adj_grphPC2))


tidygraph::as_tibble(alledges)
uniongraph
tidygraph::as_tibble(uniongraph)
uniongraph
tidygraph::as_tibble(uniongraph)
tidygraph::as_tbl_graph(uniongraph)
tidygraph::as_tbl_graph(uniongraph) |> as_tibble()
resultgraph <- graph.adjacency(get.adjacency(uniongraph))

