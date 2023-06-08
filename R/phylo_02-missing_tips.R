# match tree and meansmasses data
library(ape)
library(readr)
library(dplyr)
library(stringr)

# inputs 
rodtrees <- read.nexus("2021_conv/phylo/rodents_treeblk.nex") %>% .uncompressTipLabel()
meansmasses <- read_csv("2021_conv/data/meansMasses265.csv") %>% 
  mutate(sp=gsub(" ","_",sp))

intree <- rodtrees[[1]]$tip.label
indata <- meansmasses$sp

#tibble(missingsp=setdiff(indata,intree)) %>% write.csv(quote = F)

# edit species names in data
# missing tip actions

## edit data
mt_ed <- read_csv("2021_conv/data/missingTipActions.csv") %>% 
  filter(action=="edit data") %>% 
  mutate(indata=missingsp,fordata=intree)  
  
meansmasses_std <- 
left_join(meansmasses,mt_ed,by=c("sp"="missingsp")) %>%
  mutate(sp=if_else(sp==indata&!is.na(fordata),fordata,sp)) %>% 
  select(sp:LRa)

write_csv(meansmasses_std,"2021_conv/data/meansMasses265_rdy.csv")

# edit tree tips
mt_et <- 
  read_csv("2021_conv/data/missingTipActions.csv") %>% 
  filter(action=="edit tree tip") %>% 
  mutate(fortree=missingsp)  


## Function for renaming tips by Thomas G. and Ben Bolker on Stack Overflow
rename.tips.phylo <- function(tree, oldnames,newnames) {
  tree$tip.label[match(oldnames,tree$tip.label)] <- newnames
  return(tree)
}

# edit tips in all trees
rodtreesrn <- lapply(rodtrees, rename.tips.phylo,oldnames=mt_et$intree,
       newnames=mt_et$fortree)

class(rodtreesrn) <- "multiPhylo"
rodtreesrn

# write as Nexus file
write.nexus(rodtreesrn,file="2021_conv/phylo/rodent_trees_rdy.nex")



