library(windex)
library(purrr)
library(dplyr)
library(fs)
library(l1ou)
library(tidyr)
library(evobiR)
library(stringr)
library(tidygraph)
library(ggplot2)
library(readr)

# Load PCs cr convergent regime models ----
# PC1cr 
PC1_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC1_convRegs_crAll.rds")

# phylogenies from each model
alltreesPC1 <- map(PC1_conv_regs_cr,"tree")
# all tip labels (original order)
tips_vec <- alltreesPC1[[34]]$tip.label

# Load PC1 cr convergent regimes, Processed----
# PC1cr 
# relevant RDS files
PC1allcr <- read.csv("2021_conv/out/processed_edges/crPC1resEdges.csv") %>% select(-1)

# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n > quantile(PC1allcr$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
  rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
regMembership$focal <- 1
regMembership
}

PC1regs <- to_regComponents(PC1allcr)
split_regsPC1 <- PC1regs %>% group_split(regime)

# get dataset from model objects (any single one works)
PC1modY <- PC1_conv_regs_cr[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsPC1,right_join,PC1modY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesPC1,~ReorderData(.x,w_dfs[[1]]))
w_dfsReordered2 <-map(alltreesPC1,~ReorderData(.x,w_dfs[[2]]))

# calculate indices
W_regime1 <- map2(alltreesPC1, w_dfsReordered1, function(x = alltreesPC1[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

W_regime2 <- map2(alltreesPC1, w_dfsReordered2, function(x = alltreesPC1[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()
Wvals2 <- W_regime2 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1),
          data.frame(w=Wvals2,reg=2))
wdat$PC <- "PC1cr"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))

# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/PC1crw.csv")

# Load PC2 cr convergent regimes, Processed----
# PC2cr 
# relevant RDS files
PC2allcr <- read.csv("2021_conv/out/processed_edges/crPC2resEdges.csv") %>% select(-1)
# Load PCs cr convergent regime models ---
# PC2cr 
PC2_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC2_convRegs_crAll.rds")

# phylogenies from each model
alltreesPC2 <- map(PC2_conv_regs_cr,"tree")
# all tip labels (original order)
tips_vec <- alltreesPC2[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n > quantile(PC2allcr$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

PC2regs <- to_regComponents(PC2allcr)
split_regsPC2 <- PC2regs %>% group_split(regime)

# get dataset from model objects (any single one works)
PC2modY <- PC2_conv_regs_cr[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsPC2,right_join,PC2modY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesPC2,~ReorderData(.x,w_dfs[[1]]))
w_dfsReordered2 <-map(alltreesPC2,~ReorderData(.x,w_dfs[[2]]))

# calculate indices
W_regime1 <- map2(alltreesPC2, w_dfsReordered1, function(x = alltreesPC2[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

W_regime2 <- map2(alltreesPC2, w_dfsReordered2, function(x = alltreesPC2[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()
Wvals2 <- W_regime2 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1),
                  data.frame(w=Wvals2,reg=2))
wdat$PC <- "PC2cr"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))
wdat
# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/PC2crw.csv")

# Load PC3 cr convergent regimes, Processed----
# PC3cr 
# relevant RDS files
PC3allcr <- read.csv("2021_conv/out/processed_edges/crPC3resEdges.csv") %>% select(-1)
# Load PCs cr convergent regime models 
PC3_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/PC3_convRegs_crAll.rds")

# phylogenies from each model
alltreesPC3 <- map(PC3_conv_regs_cr,"tree")
# all tip labels (original order)
tips_vec <- alltreesPC3[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n > quantile(PC3allcr$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

PC3regs <- to_regComponents(PC3allcr)
split_regsPC3 <- PC3regs %>% group_split(regime)

# get dataset from model objects (any single one works)
PC3modY <- PC3_conv_regs_cr[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsPC3,right_join,PC3modY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesPC3,~ReorderData(.x,w_dfs[[1]]))
w_dfsReordered2 <-map(alltreesPC3,~ReorderData(.x,w_dfs[[2]]))

# calculate indices
W_regime1 <- map2(alltreesPC3, w_dfsReordered1, function(x = alltreesPC3[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

W_regime2 <- map2(alltreesPC3, w_dfsReordered2, function(x = alltreesPC3[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()
Wvals2 <- W_regime2 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1),
                  data.frame(w=Wvals2,reg=2))
wdat$PC <- "PC3cr"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))
wdat
# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/PC3crw.csv")

# Load cr multivariate convergent regimes, Processed----
# Multcr 
# relevant RDS files
Multallcr <- read.csv("2021_conv/out/processed_edges/crMultresEdges.csv") %>% select(-1)
# Load PCs cr convergent regime models 
Mult_conv_regs_cr <- readRDS("2021_conv/out/processed_l1ouModels/mult_convRegs_crAll.rds")

# phylogenies from each model
alltreesMult <- map(Mult_conv_regs_cr,"tree")
# all tip labels (original order)
tips_vec <- alltreesMult[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n > quantile(Multallcr$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

Multregs <- to_regComponents(Multallcr)
split_regsMult <- Multregs %>% group_split(regime)

# get dataset from model objects (any single one works)
MultmodY <- Mult_conv_regs_cr[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsMult,right_join,MultmodY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesMult,~ReorderData(.x,w_dfs[[1]]))

w_dfsReordered1

# calculate indices
W_regime1 <- map2(alltreesMult, w_dfsReordered1, function(x = alltreesMult[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = c("Dim.1","Dim.2","Dim.3"), focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- data.frame(w=Wvals1,reg=1)
wdat$PC <- "Multcr"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))
wdat
# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/Multcrw.csv")

############# External
# Load PC1 ext convergent regimes, Processed----
# PC1ext 
# relevant RDS files
PC1allext <- read.csv("2021_conv/out/processed_edges/extPC1resEdges.csv") %>% select(-1)
# Load PCs ext convergent regime models 
PC1_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC1_convRegs_extAll.rds")

# phylogenies from each model
alltreesPC1 <- map(PC1_conv_regs_ext,"tree")
# all tip labels (original order)
tips_vec <- alltreesPC1[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n > quantile(PC1allext$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

PC1regs <- to_regComponents(PC1allext)
split_regsPC1 <- PC1regs %>% group_split(regime)

# get dataset from model objects (any single one works)
PC1modY <- PC1_conv_regs_ext[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsPC1,right_join,PC1modY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesPC1,~ReorderData(.x,w_dfs[[1]]))
w_dfsReordered2 <-map(alltreesPC1,~ReorderData(.x,w_dfs[[2]]))
w_dfsReordered3 <-map(alltreesPC1,~ReorderData(.x,w_dfs[[3]]))

# calculate indices
W_regime1 <- map2(alltreesPC1, w_dfsReordered1, function(x = alltreesPC1[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

W_regime2 <- map2(alltreesPC1, w_dfsReordered2, function(x = alltreesPC1[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

W_regime3 <- map2(alltreesPC1, w_dfsReordered3, function(x = alltreesPC1[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})
# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()
Wvals2 <- W_regime2 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()
Wvals3 <- W_regime3 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1),
                  data.frame(w=Wvals2,reg=2),
                  data.frame(w=Wvals3,reg=3))
wdat$PC <- "PC1ext"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))

# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/PC1extw.csv")

# Load PC2 ext convergent regimes, Processed----
# PC2ext 
# relevant RDS files
PC2allext <- read.csv("2021_conv/out/processed_edges/extPC2resEdges.csv") %>% select(-1)
# Load PCs ext convergent regime models 
PC2_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC2_convRegs_extAll.rds")

# phylogenies from each model
alltreesPC2 <- map(PC2_conv_regs_ext,"tree")
# all tip labels (original order)
tips_vec <- alltreesPC2[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n >= quantile(PC2allext$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

PC2regs <- to_regComponents(PC2allext)
split_regsPC2 <- PC2regs %>% group_split(regime)

# get dataset from model objects (any single one works)
PC2modY <- PC2_conv_regs_ext[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsPC2,right_join,PC2modY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesPC2,~ReorderData(.x,w_dfs[[1]]))

# calculate indices
W_regime1 <- map2(alltreesPC2, w_dfsReordered1, function(x = alltreesPC2[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1))
wdat$PC <- "PC2ext"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))

# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/PC2extw.csv")

# Load PC3 ext convergent regimes, Processed----
# PC3ext 
# relevant RDS files
PC3allext <- read.csv("2021_conv/out/processed_edges/extPC3resEdges.csv") %>% select(-1)
# Load PCs ext convergent regime models 
PC3_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC3_convRegs_extAll.rds")

# phylogenies from each model
alltreesPC3 <- map(PC3_conv_regs_ext,"tree")
# all tip labels (original order)
tips_vec <- alltreesPC3[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n >= quantile(PC3allext$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

PC3regs <- to_regComponents(PC3allext)
split_regsPC3 <- PC3regs %>% group_split(regime)

# get dataset from model objects (any single one works)
PC3modY <- PC3_conv_regs_ext[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsPC3,right_join,PC3modY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesPC3,~ReorderData(.x,w_dfs[[1]]))

# calculate indices
W_regime1 <- map2(alltreesPC3, w_dfsReordered1, function(x = alltreesPC3[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = "V1", focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1))
wdat$PC <- "PC3ext"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))

# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/PC3extw.csv")

# Load ext multivariate convergent regimes, Processed----
# Multext 
# relevant RDS files
Multallext <- read.csv("2021_conv/out/processed_edges/extMultresEdges.csv") %>% select(-1)
# Load PCs cr convergent regime models 
Mult_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/mult_convRegs_extAll.rds")

# phylogenies from each model
alltreesMult <- map(Mult_conv_regs_ext,"tree")
# all tip labels (original order)
tips_vec <- alltreesMult[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n > quantile(Multallext$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

Multregs <- to_regComponents(Multallext)
split_regsMult <- Multregs %>% group_split(regime)

# get dataset from model objects (any single one works)
MultmodY <- Mult_conv_regs_ext[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsMult,right_join,MultmodY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesMult,~ReorderData(.x,w_dfs[[1]]))
w_dfsReordered2 <-map(alltreesMult,~ReorderData(.x,w_dfs[[2]]))

# calculate indices
W_regime1 <- map2(alltreesMult, w_dfsReordered1, function(x = alltreesMult[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = c("Dim.1","Dim.2","Dim.3"), focal = .y$focal))
})
W_regime2 <- map2(alltreesMult, w_dfsReordered2, function(x = alltreesMult[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = c("Dim.1","Dim.2","Dim.3"), focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()
Wvals2 <- W_regime2 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- bind_rows(data.frame(w=Wvals1,reg=1),
                  data.frame(w=Wvals2,reg=2))
wdat$PC <- "Multext"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))
wdat
# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/Multextw.csv")

# Load ALLTRaits multivariate convergent regimes, Processed----
# MultALL (wb for whole body)
# relevant RDS files
Multallwb <- read.csv("2021_conv/out/processed_edges/allPCsMultresEdges.csv") %>% select(-1)
# Load PCs cr convergent regime models 
Mult_conv_regs_all <- readRDS("2021_conv/out/processed_l1ouModels/ALLPCs_convRegs_All.rds")

# phylogenies from each model
alltreesMultwb <- map(Mult_conv_regs_all,"tree")
# all tip labels (original order)
tips_vec <- alltreesMultwb[[34]]$tip.label


# split into components, one per regime
to_regComponents <- function(PCsummary) {
  PCsummary %>%
    filter(n >= quantile(Multallwb$n,0.90)) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$focal <- 1
  regMembership
}

Multregswb <- to_regComponents(Multallwb)
split_regsMultwb <- Multregswb %>% group_split(regime)

# get dataset from model objects (any single one works)
MultwbmodY <- Mult_conv_regs_all[[34]]$Y %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sp")

# join Y with conv regimes to tag convergent taxa
# for windex
w_dfs <- map(split_regsMultwb,right_join,MultwbmodY)

# set up logical vector to indicate convergent and non-conv taxa
w_dfs <- w_dfs %>% map(replace_na,list(focal=0))
# rename for windex
w_dfs <- w_dfs %>% map(~rename(.x,species=sp))

# add rownames for reordering
dfsetrownames <- function(targetdf,which_col){
  targetdf <- distinct(targetdf)
  rownames(targetdf) <- targetdf[[which_col]]
  targetdf
}
w_dfs <- w_dfs %>% map(dfsetrownames,"species")

# reorder data to match tree tips
length(w_dfs)

w_dfsReordered1 <-map(alltreesMultwb,~ReorderData(.x,w_dfs[[1]]))

# calculate indices
W_regime1 <- map2(alltreesMultwb, w_dfsReordered1, function(x = alltreesMultwb[1], y = w_dfsReordered1[1]) {
  map2(list(x), list(y), ~windex(data.frame(.y),.x,traits = c("Dim.1","Dim.2","Dim.3","Dim.4","Dim.5"), focal = .y$focal))
})

# extract w as vectors
Wvals1 <- W_regime1 %>% map_depth("w",.depth = 2) %>% flatten() %>% flatten_dbl()

wdat <- data.frame(w=Wvals1,reg=1)
wdat$PC <- "MultALLPCs"

ggplot(wdat)+
  geom_boxplot(aes(factor(reg),w))
wdat
# export
write_csv(wdat,file = "2021_conv/out/processed_regWs/MultALLPCsw.csv")


#### junk below this

W_regime2 %>% map_depth("w",.depth = 2) 
map2(w_dfsReordered1,alltreesPC1, ~windex(.y,.x,traits = "V1", focal = .y$focal))


windex(data.frame(w_dfsReordered1[[1]]),alltreesPC1[[1]],traits = "V1", 
       focal =w_dfsReordered1[[1]]$focal) 
w_dfsReordered1[[1]]


ww <- map_depth(W_regime1,.depth = 2,"w") %>% purrr::simplify() %>% 
  map(`[`,1)

ww %>% flatten() %>% flatten()
W_regime1 %>% tibble::enframe() %>% unnest_auto(value)
length(alltreesPC1)
length(w_dfsReordered1)

w_dfs
# # hack for map2
# w_dfs <- w_dfs %>% map(replicate,n=length(alltreesPC1),simplify=FALSE)
# w_dfs <- w_dfs %>% purrr::simplify()
# w_dfs
# reorder data to match tree tips
length(w_dfs)

w_dfsReordered <- list()
w_dfs[1]

for(i in seq_along(1:length(w_dfs))){
w_dfsReordered[i] <- map(alltreesPC1,~ReorderData(.x,w_dfs[[i]]))  
}

w_dfsReordered
map(alltreesPC1[1:3],~ReorderData(.x,w_dfs[[1]]))
map(alltreesPC1[1:3],~ReorderData(.x,w_dfs[[2]]))

alltreesPC1
w_dfs[2]


map2(alltreesPC1, w_dfsReordered1, function(x = alltreesPC1[1], y = w_dfsReordered1[[1]]) {
map2(x, y, ~windex(.y,.x,traits = "V1", focal = .y$focal))
     })

# convregDFs <- 
map2(alltreesPC1,w_dfs, 
     function(x = alltreesPC1[[1]], y = w_dfs[[1]][[1]]) {
    map2(list(x), y, ReorderData)
})

imap(alltreesPC1,ReorderData,w_dfs)

map2(alltreesPC1,w_dfs[1],ReorderData)
