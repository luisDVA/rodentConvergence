# graft missing species to tree
library(ape)
library(readr)
library(dplyr)
library(treeman)
library(stringr)
library(purrr)

# inputs 
rodtrees <- read.nexus("2021_conv/phylo/rodent_trees_rdy.nex")
# edit tree tips
mt_et <- 
  read_csv("2021_conv/data/missingTipActions.csv") %>% 
  filter(stringr::str_detect(action,"graft")) %>% 
  mutate(fortree=missingsp)  
mt_et


# define fns

# edge lengths for all congeners across all trees
getGenusAges <- function(treelist, regpattern) {
  # from phytools blog
  gentips <- treelist$tip.label[which(str_detect(treelist$tip.label, regpattern))]
  gennodes <- sapply(gentips, function(x, y) which(y == x), y = treelist$tip.label)
  treelist$edge.length[sapply(gennodes, function(x, y) which(y == x), y = treelist$edge[, 2])]
}

get_tip_age <- function(tree,tipname){
  tmtree <- as(tree,'TreeMan')
  getSpnAge(tmtree,tipname,getAge(tmtree))
}  


add_tips_multi <- function(tree,newtip,targettip,agesvec){
  tmtree <- as(tree,'TreeMan')  
  tmage <-  getAge(tmtree)
  target.tip.age <- getSpnAge(tmtree,targettip,getAge(tmtree))$start
  new.tip.age <- 100
  while (new.tip.age > target.tip.age){
    new.tip.age <- sample(agesvec,1)
  }
  newtreetm <- addTip(tmtree,tid = newtip,sid = targettip,
         strt_age = new.tip.age,end_age = 0,tree_age = tmage)
  newtree <- as(newtreetm, 'phylo')  
  newtree$node.label <- NULL
  newtree
  }

# Abrothrix hirta grafted to A. longipilis
abrothrix_edges <- rodtrees %>% map(getGenusAges,"(?!Abrothrix_longipilis)Abrothrix") %>%
  flatten_dbl()

rodg1 <-  
rodtrees %>% 
  map(~add_tips_multi(.x,newtip = "Abrothrix_hirta",
                                        targettip = "Abrothrix_longipilis",
                                        agesvec = abrothrix_edges))

# Octodon ricardojeda grafted to O. bridgesi
octodon_edges <- rodg1 %>% map(getGenusAges,"(?!Octodon_bridgesi)Octodon") %>%
  flatten_dbl()
rodg2 <-  
  rodg1 %>% 
  map(~add_tips_multi(.x,newtip = "Octodon_ricardojeda",
                      targettip = "Octodon_bridgesi",
                      agesvec = octodon_edges))

# Lagidium moreni grafted to L. viscacia
lagidium_edges <- rodg2 %>% map(getGenusAges,"(?!Lagidium_viscacia)Lagidium") %>%
  flatten_dbl()

rodtreesALL <-  
  rodg2 %>% 
  map(~add_tips_multi(.x,newtip = "Lagidium_moreni",
                      targettip = "Lagidium_viscacia",
                      agesvec = lagidium_edges))

rodtreesALL
class(rodtreesALL) <- "multiPhylo"

rodtreesALL
# write as Nexus file
write.nexus(rodtreesALL,file="2021_conv/phylo/rodents_treeblk_rdy.nex")


