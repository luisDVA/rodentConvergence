# match tree and meansmasses data
library(ape)
library(readr)
library(purrr)

# inputs 
rodtrees <- read.nexus("2021_conv/phylo/rodents_treeblk_rdy.nex")
class(rodtrees)
meansmasses <- read_csv("2021_conv/data/meansMasses265_rdy.csv") 

# subset to data
phydata <- lapply(rodtrees,keep.tip,tip=meansmasses$sp) 
class(phydata) 
class(phydata) <- "multiPhylo"
phydata[[1]]
plot(phydata[[1]])
is.ultrametric(phydata[[1]])
phydata <- map(phydata,~phytools::force.ultrametric(.x,message=FALSE ))
class(phydata) <-  'multiPhylo'
# export phylo for l1ou
write.nexus(phydata,file="2021_conv/phylo/rodent_trees265.nex")
