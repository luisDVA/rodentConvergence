#############################################################################
#summarize.lasso.convergence.output FUNCTION FROM SLATER & FRISCIA (2019)
summarize.lasso.convergence.output <- function(model) {
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
}

regimes_summary <- function(model){
  n_shifts <- model$nShifts
  sc <- model$shift.configuration
  regs <- sort(unique(names(sc)))
  regimesls <- map(regs,~sc[which(names(sc) == .x)])
  n_regs <- length(regimesls)
  n_conv <- sum(lengths(regimesls)>1)
  tibble(shifts=n_shifts,regimes=n_regs,conv_regimes=n_conv)
}
