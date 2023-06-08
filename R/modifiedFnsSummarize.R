#############################################################################
# modified FUNCTION FROM SLATER & FRISCIA (2019)

summarize.lasso.convergence.output <- function(model) {
  # check for convergence first
  n_shifts <- model$nShifts
  sc <- model$shift.configuration
  regs <- sort(unique(names(sc)))
  regimesls <- map(regs,~sc[which(names(sc) == .x)])
  n_regs <- length(regimesls)
  n_conv <- sum(lengths(regimesls)>1)
  if (n_conv==0) return(data.frame(focal=c(),sp=c(),n=c()))
  
  
  bt <- branching.times(model$tree)
  shift.times.max <- bt[match(model$tree$edge[model$shift.configuration, 1], names(bt))]
  shift.crown <- model$tree$edge[model$shift.configuration, 2]
  
  shift.times.min <- bt[match( shift.crown, names(bt))]
  
  shift.taxa <- list()
  
  for(i in 1:length(shift.crown)) shift.taxa[[i]]<- (tips(model$tree, shift.crown[i]))
  
  ogshift.taxa <- shift.taxa
  taxa_regimes <- map_df(ogshift.taxa,as.data.frame,.id = 'regime')
  names(taxa_regimes) <- c("og_reg","sp")
  
  names(shift.taxa) <- names(model$shift.configuration)
  ordered.list <- list()
  list.order <- order(as.numeric(unlist(lapply(shift.taxa, length))))
  for(i in 1:length(list.order)) {
    ordered.list[[i]] <- shift.taxa[[list.order[i]]]
  }
  
  names(ordered.list) <- names(shift.taxa)[list.order]
  
  for(i in 1:(length(ordered.list)-1)) {
    for(j in (i+1): (length(ordered.list))) {
      if(!all(is.na(match(ordered.list[[i]], ordered.list[[j]])))) {
        ordered.list[[j]] <- ordered.list[[j]][-match(ordered.list[[i]], ordered.list[[j]])]
      }
    }
  }
  
  convergent.regimes <- list()
  new.regimes<-unique(names(ordered.list))
  for(i in 1:length(new.regimes)){
    xx <- which(names(ordered.list)==new.regimes[[i]])
    tmp <- character()
    for(j in 1:length(xx)) {
      tmp <-  append(tmp, ordered.list[[xx[j]]])
    }
    convergent.regimes[[i]] <-tmp
    rm(tmp)
  }
  
  convergent.regimes[[length(new.regimes)+1]] <- setdiff(model$tree$tip.label, unlist(convergent.regimes)) 
  names(convergent.regimes) <- c(new.regimes, "noshift")
  convergent.regimes  
  as.character(unlist(convergent.regimes))
  reg <- names(convergent.regimes)
  yy <- character()
  for(i in 1:length(reg)) {
    yy<-append(yy, rep(reg[i], length(convergent.regimes[[i]])))
  }
  
  convregdf <- as.data.frame(cbind(as.character(unlist(convergent.regimes)), yy))
  convregdf
  names(convregdf) <- c("sp","conv_reg")
  convregdf <- left_join(convregdf,taxa_regimes)
  convregdf$og_reg <- ifelse(is.na(convregdf$og_reg),"noshift",convregdf$og_reg)
  convregdf
  #tag 
  sc <- model$shift.configuration
  regs <- sort(unique(names(sc)))
  regimesls <- map(regs,~sc[which(names(sc) == .x)])
  n_regs <- length(regimesls)
  c_regs <- regs[lengths(regimesls)>1]
  cregsdf <- data.frame(reg=c_regs)
  cregsdf$regime_type="conv_reg"
  convregdf <- left_join(convregdf,cregsdf,by=c("conv_reg"="reg"))
  convregdf
  convregdf$regime_type=ifelse(is.na(convregdf$regime_type),"nc",convregdf$regime_type)
  convregdf  
  }


# count shifts, regimes, and convergent regimes
regimes_summary <- function(model){
  n_shifts <- model$nShifts
  sc <- model$shift.configuration
  regs <- sort(unique(names(sc)))
  regimesls <- map(regs,~sc[which(names(sc) == .x)])
  n_regs <- length(regimesls)
  n_conv <- sum(lengths(regimesls)>1)
  tibble(shifts=n_shifts,regimes=n_regs,conv_regimes=n_conv)
}

# summarize convergent regimes across a list of models
summarise_convregs <- function(convRegList) {
  # condense convergent clades
  convtipsPC <- map(convRegList, summarize.lasso.convergence.output)
  convtipsPC <- discard(convtipsPC, ~nrow(.) == 0)
  convtipslongPC <- convtipsPC %>% map_df(bind_rows, .id = "tree")
  justcon_PC <- convtipsPC %>% map(~ filter(.x, conv_reg != "noshift"))
  justcon_PClong <- justcon_PC %>% map_df(bind_rows, .id = "tree")
  
  
  # all species
  allsps_conv <- convtipslongPC %>%
    distinct(sp) %>%
    pull(sp)
  
  # 'basal' species across all trees
  treesgp <- convtipslongPC %>%
    group_by(tree) %>%
    n_groups()
  noShift_taxa <- convtipslongPC %>%
    group_by(tree) %>%
    filter(conv_reg == "noshift") %>%
    ungroup() %>%
    count(sp) %>%
    filter(n == treesgp) %>%
    pull(sp)
  
  convtaxa <- setdiff(allsps_conv, noShift_taxa)
  
  
  # get all convergent taxa for a focal sp, remove sps from same regime
  get_conv_taxa <- function(modeldf, species) {
    mbros <- modeldf %>% filter(stringr::str_detect(sp, species))
    n_basal <- nrow(mbros[mbros$conv_reg == "noshift", ])
    mbros_convOnly <- mbros %>% filter(conv_reg != "noshift")
    conv_regimeN <- pull(mbros_convOnly, conv_reg)
    og_regimeN <- pull(mbros_convOnly, og_reg)
    convBros <- modeldf %>% filter(conv_reg %in% conv_regimeN)
    conv_sps <-
      convBros %>%
      filter(sp != species) %>%
      filter(!og_reg %in% og_regimeN) %>%
      mutate(focal = species, .before = 1) %>%
      select(focal, sp) %>%
      distinct()
    return(list(noshift = n_basal, convergent_sps = conv_sps))
  }
  
  #
  count_bros <- function(focal_sp) {
    convtipsPC %>%
      map(~ get_conv_taxa(.x, focal_sp)) %>%
      map_df("convergent_sps") %>%
      add_count(sp)
  }
  
  PCallbros <- map(convtaxa, count_bros) %>% map_df(bind_rows)
  PCallbros <- PCallbros %>% distinct()
  PCallbros
}