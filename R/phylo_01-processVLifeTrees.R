# # trim up2019 (v4) phylo to meansmasses265 dataset
library(ape)
library(fs)
library(evobiR)
library(stringr)
library(dplyr)


# sample 100 random trees 
SampleTrees(trees = "/media/luisd/H/vertlife_mammal_phylo/doi_10.5061_dryad.tb03d03__v4/Data_S7_Mammalia_credibleTreeSets_tipDR/Data_S7_Mammalia_credibleTreeSets_tipDR/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_nexus.trees",
                            burnin=0.2,format="nex",final.number = 100,prefix="samp_")

# move to correct path
fs::file_move("samp_ .nex",new_path = "2021_conv/phylo/")
fs::file_move("2021_conv/phylo/samp_ .nex","2021_conv/phylo/tree_block.nex")

# all mammals
treeblock <- read.nexus("2021_conv/phylo/tree_block.nex")

# susbset order
alltips <- treeblock[[1]]$tip.label

# rodents
rodent_tipIndices <- which(str_detect(alltips,"RODENTIA$"))
rodent_tips <- alltips[rodent_tipIndices]

# trim tree block
# define keep.tip function by Liam Revell
keep.tip<-function(tree,tip) drop.tip(tree,setdiff(tree$tip.label,tip))

# trim
RodTreesF <-lapply(treeblock,keep.tip,tip=rodent_tips)
class(RodTreesF) <- "multiPhylo"

drophigher <- function(tree){
tree$tip.label <- str_extract(tree$tip.label,".*(?=(_.*?_.*))")  
tree
}

RodTrees_rdy <- lapply(RodTreesF,drophigher)


class(RodTrees_rdy) <- "multiPhylo"

# write as Nexus file
write.nexus(RodTrees_rdy,file="2021_conv/phylo/rodents_treeblk.nex")





