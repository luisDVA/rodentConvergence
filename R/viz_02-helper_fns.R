# edges plots helpers
# functions to prepare and plot sets of convergent species
## as networks and on phylogenetic trees

to_regComponents <- function(PCsummary,comparison,wPC,wmtype) {
  fun1 <- match.fun(comparison)
  PCsummary %>%
    filter(fun1(n,quantile(PCsummary$n,0.90))) %>%
    rename(to = 1, from = 2, weight = 3) %>%
    mutate(weight = weight) %>%
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) |> 
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$PC <- wPC
  regMembership$mtype <- wmtype
  list(regMembership,adj_grphPC)
}

# Special case for whole body data
to_regComponentswb <- function(PCsummary,wPC,wmtype) {
  PCsummary %>%
    rename(to = 1, from = 2) %>% 
    mutate(across(where(is.character),~str_replace(.x,"_"," "))) %>%  
    as_tbl_graph(directed=FALSE) -> adj_grphPC
  regMembership <- components(adj_grphPC)$membership |> tibble::enframe() |> 
    rename(sp=1,regime=2) |> mutate(sp=str_replace(sp," ","_"))
  regMembership$PC <- wPC
  regMembership$mtype <- wmtype
  list(regMembership,adj_grphPC)
}

edgesToNodes <- function(PCgrph) {
  # Separate out edges and node data frames
  tg_nodes <-
    PCgrph %>%
    activate(nodes) %>%
    data.frame() %>%
    tibble::rownames_to_column("rowid") %>%
    mutate(rowid = as.integer(rowid))
  tg_edges <-
    PCgrph %>%
    activate(edges) %>%
    data.frame()
  named_edge_list <-
    tg_edges %>%
    # Rename from nodes
    left_join(tg_nodes, by = c("from" = "rowid")) %>%
    select(-from) %>% # Remove unneeded column
    rename(from = name) %>% # Rename column with names now
    # Rename to nodes
    left_join(tg_nodes, by = c("to" = "rowid")) %>%
    select(-to, -w.y) %>% # Remove unneeded column
    rename(to = name, w = w.x, ntrees = weight)
  named_edge_list %>% as_tbl_graph(directed = FALSE)
}
makeWgraph <- function(regsobj, PClab) {
  regsobj[[2]] %>%
    activate(nodes) %>%
    mutate(reg = group_components(), PC = PClab) %>%
    left_join(regWs_summary) %>%
    edgesToNodes()
}

makeWgraphwb <- function(regsobj, PClab) {
  regsobj[[2]] %>%
    activate(edges) %>%
    mutate(reg = treeID, PC = "MultALLPCs") %>%
    mutate(reg= as.numeric(factor(treeID))) %>% 
    left_join(regWs_summary) 
}

pltnetwork <- function(regsObjectgrp, lyt = "kk") {
  ncolorssc <- regsObjectgrp %>%
    activate(edges) %>%
    as_tibble() %>%
    pull(ntrees) %>%
    n_distinct()
  ggraph(regsObjectgrp, layout = lyt) +
    geom_edge_link(aes(
      color = factor(ntrees),
      alpha = factor(ntrees),
      edge_width = factor(round(w, 2))
    )) +
    scale_edge_width_discrete(range = c(0.3, 1.3), name = "Wheatsheaf Index") +
    scale_edge_alpha_discrete(
      range = c(0.3, 0.6),
      name = "Number of\n trees"
    ) +
    geom_node_point(size = 0.7,pch=21) +
    scale_edge_colour_manual(
      values = scico(ncolorssc,
        palette = "devon", begin = 0.1, end = 0.7,
        direction = -1
      ), name = "Number of\n trees"
    ) +
    geom_text_repel(aes(x, y, label = name,segment.inflect = TRUE), family="Laksaman",
                    segment.size=0.2,segment.color="#4a4e69",
                    fontface = "italic", color = "black", bg.color = "white", bg.r = 0.15,size=2.2) +
    theme_graph()+theme(legend.position = "bottom",
                        legend.text = element_text(family = "Laksaman"))
}

pltnetworkwb <- function(regsObjectgrp, lyt = "kk") {
  pcolors <- c("#001870FF","#F85820FF","#686868FF","#006ade","#ff8b83","#62bd41")
  ggraph(regsObjectgrp, layout = lyt) +
    geom_edge_link(aes(
      color = factor(reg),group=factor(reg),
      edge_width = factor(round(w, 2))
    )) +
    scale_edge_width_discrete(range = c(0.3, 1.3), name = "Wheatsheaf Index") +
    geom_node_point(size = 0.7, pch=21) +
    scale_edge_colour_manual(
      values = pcolors, 
      name = "Tree"
    ) +
    geom_text_repel(aes(x, y, label = name,segment.inflect = TRUE), family="Laksaman",
                    segment.size=0.2,segment.color= "#4a4e69",
                    fontface = "italic", color = "black", bg.color = "white", bg.r = 0.15,size=2.2) +
    theme_graph()+theme(legend.position = "bottom",
                        legend.text = element_text(family = "Laksaman"))
}

# phylo steps
# get otu data and join to tree
join_groupINFO <- function(regssObj) {
  groupsInfo <- split(regssObj[[1]]$sp, regssObj[[1]]$regime)
  # join Regs with the tree
  groupOTU(rodphyl, groupsInfo, overlap = "abandon")
}

# for wb model, use initial edges table
join_groupINFOwb <- function(edgesDat) {
  Multregslong <- edgesDat %>%
    pivot_longer(-treeID) %>% 
    rename(sp=value) %>% select(-name) %>% 
    distinct(sp,treeID) %>% 
    select(sp,treeID) 
  regssObj <- Multregslong 
  groupsInfo <- 
    split(regssObj$sp, regssObj$treeID)
  # join Regs with the tree
  groupOTU(rodphyl, groupsInfo, overlap = "abandon")
}

# process groups and plot phylogeny
convPhyloPlot <- function(groupedPhyl) {
  # tidytree manipulation steps
  otutib <- groupedPhyl %>% as_tibble()
  otutibsubset <- otutib %>%
    filter(group != "0") %>%
    filter(!is.na(label))
  # get higher taxa
  highertaxPhyl <- left_join(otutibsubset, spMDDnamesTax, by = c("label" = "sp"))
  highertaxPhyl <-
    highertaxPhyl %>%
    mutate(genustiplab = str_extract(label, "^[^_]+")) %>%
    add_count(genustiplab, name = "tiplab_congeners") %>%
    add_count(tribe, name = "tribemembers") %>%
    mutate(
      tribemembers =
        case_when(
          is.na(tribe) ~ NA_integer_,
          TRUE ~ tribemembers
        )
    ) %>%
    add_count(family, name = "fam_members")
  # for species level labels
  highertaxPhyl <-
    highertaxPhyl %>%
    mutate(short_label = str_replace(label, "_", " ")) %>%
    tidyr::separate(short_label, into = c("gen", "spep")) %>%
    mutate(gen = str_replace(gen, "\\B[a-z]+", "\\.")) %>%
    unite("short_label", gen:spep, sep = " ")

  # final clade labels
  highertaxPhyl <- highertaxPhyl %>%
    mutate(cl_label = coalesce(tribe, subfamily, family)) %>%
    mutate(cl_label = ifelse(fam_members == 1, short_label, cl_label)) %>%
    mutate(cl_label = str_to_sentence(cl_label))

  # process nodes for label placement
  formrca <- highertaxPhyl %>%
    select(cl_label, node) %>%
    add_count(cl_label, name = "numnodes")
  formrca <- formrca %>% filter(numnodes != 1)
  singletons <- highertaxPhyl %>%
    add_count(cl_label, name = "numnodes") %>%
    filter(numnodes == 1) %>%
    select(label, cl_label) %>%
    mutate(label2 = cl_label)
  nodens <- formrca %>%
    group_split(cl_label) %>%
    purrr::map(pull, node)
  gnamesvec <- formrca %>%
    group_split(cl_label) %>%
    purrr::map(pull, cl_label) %>%
    purrr::map(unique)
  groupedPhyl <- left_join(groupedPhyl, singletons)
  # define colors once
  groupcolors <- c("#bababa", "#e7298a", "#315a7b",
                    "#ff5a5a","#b44141","#414a83",
                    "#7b7bb4")

  pcfig <-
    ggtree(groupedPhyl, layout = "daylight", aes(color = group)) +
    scale_color_manual(
      values = groupcolors,
      guide = "none"
    ) +
    theme(legend.position = "right")

  if (nrow(formrca)==0) {
    pcfig <- pcfig +
      geom_tiplab(aes(label=label2,color=group),
                  family="Laksaman",size=2.3,hjust = -0.1)
    pcfig
  }
  
  else {
  
  # define function for getting MRCA node
  getCladeNode <- function(ftreej, nodesvec, gname) {
    nodenum <- getMRCA(ftreej, tip = nodesvec)
    tibble(clade = gname, node = nodenum)
  }

  # table with genera and node numbers
  genNodes <-
    purrr::map2_df(nodens, gnamesvec, ~ getCladeNode(as.phylo(groupedPhyl), .x, .y)) # %>%
  nodeslabs <- bind_rows(genNodes, singletons)
  # loop over colors too
  genNodes <- genNodes %>%
    left_join(as_tibble(groupedPhyl)) %>%
    select(node, clade, group)
  colorsdf <- data.frame(
    group = levels(genNodes$group),
    hexcol = groupcolors[1:length(levels(genNodes$group))]
  )
  genNodes <- genNodes %>% left_join(colorsdf)

  # loop to draw the gen labels
  for (j in 1:nrow(genNodes)) {
    # Then add each clade label
    pcfig <- pcfig +
      geom_cladelab(
        node = genNodes$node[j],
        label = genNodes$clade[j],
        textcolor = genNodes$hexcol[j],
        barsize = 0.1, offset = -0.1,
        offset.text = 1,
        vjust = 0.5, angle = "auto",
        fontsize = 2.3,family="Laksaman"
      )
  }

  # add non-clade nodes
  pcfig <- pcfig +
    geom_tiplab(aes(label = label2, color = group), 
                size = 2.3, hjust = -0.1,family="Laksaman")

  pcfig
  }
}
