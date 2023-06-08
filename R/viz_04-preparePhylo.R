# draw rodent phylogeny, needs objects from previous viz scripts
rodphyl
groupedPhyl <- rodphyl
# process groups and plot phylogeny
# tidytree manipulation steps
  otutib <- groupedPhyl %>% as_tibble()
  otutib
  otutibsubset <- otutib

# get higher taxa
highertaxPhyl <- left_join(otutibsubset, spMDDnamesTax, by = c("label" = "sp"))
  
# final clade labels
  highertaxPhyl <- highertaxPhyl %>%
    mutate(cl_label = coalesce(subfamily, family)) %>%
    mutate(cl_label = str_to_sentence(cl_label)) %>% 
    mutate(cl_label=str_replace(cl_label,"Incertae sedis","Sigmodontinae"))

# change to families for some clades
highertaxPhyl <- 
  highertaxPhyl %>% 
  mutate(cl_label=if_else(family=="HETEROMYIDAE","Heteromyidae",cl_label,missing=cl_label)) %>% 
  mutate(cl_label=if_else(superfamily=="OCTODONTOIDEA","Octodontoidea",cl_label,missing=cl_label))
  
# for direct label
  fordlab <- highertaxPhyl %>%
    select(label,cl_label, node) %>%
    add_count(cl_label, name = "numnodes") %>% 
    filter(numnodes==1)
    
# process nodes for label placement
  formrca <- highertaxPhyl %>%
    select(cl_label, node) %>%
    add_count(cl_label, name = "numnodes") %>% 
    filter(!is.na(cl_label))

  
  
nodens <- formrca %>%
    group_split(cl_label) %>%
    purrr::map(pull, node)
gnamesvec <- formrca %>%
    group_split(cl_label) %>%
    purrr::map(pull, cl_label) %>%
    purrr::map(unique)

# for space, xerinae 
gnamesvec
nodens
nodens <- 
  nodens[-which(gnamesvec=="Xerinae")]
gnamesvec <-   gnamesvec[gnamesvec!="Xerinae"]

# define function for getting MRCA node
    getCladeNode <- function(ftreej, nodesvec, gname) {
      nodenum <- getMRCA(ftreej, tip = nodesvec)
      tibble(clade = gname, node = nodenum)
    }

# table with genera and node numbers
genNodes <-
    purrr::map2_df(nodens, gnamesvec, ~ getCladeNode(as.phylo(groupedPhyl), .x, .y)) 


# loop over colors too
    genNodes <- genNodes %>%
      left_join(as_tibble(groupedPhyl)) %>%
      select(node, clade) %>% na.omit()

pcfig <- ggtree(groupedPhyl,layout = "roundrect",size=0.2)
pcfig

# loop to draw the labels
  for (j in 1:nrow(genNodes)) {
      # Then add each clade label
      pcfig <- pcfig +
        geom_cladelab(
          node = genNodes$node[j],
          label = genNodes$clade[j],
         barsize = 0.7, offset = -0.6,
          offset.text = 0.1,
          vjust = 0.5, 
          fontsize = 3,family="Laksaman"
        )
    }
pcfig  

fordlab <-  
  fordlab %>% filter(!cl_label %in% c("Zapodidae","Euchoreutinae",
                                      "Erethizontinae","Dolichinae")) 
cavlab <- fordlab %>% filter(cl_label=="Caviinae")
geomylab <- fordlab%>% filter(cl_label=="Geomyinae")
castlab <- fordlab %>% filter(cl_label=="Castoridae") 

# draw the singleton labels
pcfig <-  
pcfig+geom_cladelab(node = cavlab$node,label = cavlab$cl_label,
  barsize = 0.7, offset = -0.6,
  offset.text = 0.1,
  vjust = 0.7, 
  fontsize = 3,family="Laksaman")+
  geom_cladelab(node = geomylab$node,label = geomylab$cl_label,
                      barsize = 0.7, offset = -0.6,
                      offset.text = 0.1,
                      vjust = 0.8, 
                      fontsize = 3,family="Laksaman")+
  geom_cladelab(node = castlab$node,label = castlab$cl_label,
                barsize = 0.7, offset = -0.6,
                offset.text = 0.1,
                vjust = 0.5, 
                fontsize = 3,family="Laksaman")




pcfig <- pcfig+ hexpand(0.3)
pcfig

