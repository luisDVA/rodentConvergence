# external, all 100 trees BIC, model data
library(fs)
library(l1ou)
library(purrr)
library(geiger)
library(dplyr)
library(stringr)
source("2021_conv/R/modifiedFnsSummarize.R")
library(ggplot2)
library(ggridges)
library(tidyr)

# Load PC1 ext convergent regimes----
# releant RDS files
# Load PC1 ext convergent regimes----
PC1_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC1_convRegs_extAll.rds")
PC1extRegimes <- PC1_conv_regs_ext %>% map_df(regimes_summary) %>% mutate(mod="PC1")

# Load PC2 ext convergent regimes----
PC2_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC2_convRegs_extAll.rds")
PC2extRegimes <- PC2_conv_regs_ext %>% map_df(regimes_summary) %>% mutate(mod="PC2")

# Load PC3 ext convergent regimes----
PC3_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/PC3_convRegs_extAll.rds")
PC3extRegimes <- PC3_conv_regs_ext %>% map_df(regimes_summary)%>% mutate(mod="PC3")

# multivariate ext conv regs
allPCs_conv_regs_ext <- readRDS("2021_conv/out/processed_l1ouModels/mult_convRegs_extAll.rds")
allPCsextRegimes <- allPCs_conv_regs_ext %>% map_df(regimes_summary)%>% mutate(mod="multivariate")


regresultsEXT <- bind_rows(PC1extRegimes,PC2extRegimes,PC3extRegimes,allPCsextRegimes)

regs_summary <- 
  regresultsEXT %>% group_by(mod) %>% summarise(nshifts=median(shifts),
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
  select(mod,nshifts,nregimes,nconv_regimes) %>% gt::gt()

# all shifts 
PC1_conv_regs_cr %>% map_int('nShifts') %>% median
PC1_conv_regs_cr %>% map_int('nShifts') %>% sd %>% round(2)
PC1_conv_regs_cr %>% map_int('nShifts')  %>% min
PC1_conv_regs_cr %>% map_int('nShifts')  %>% max




regresultsCR

regresultsCR %>% pivot_longer(-mod) %>% 
  ggplot(aes(y=name,x=value,height = stat(density)))+
  geom_density_ridges(stat="binline",bins=8,scale=0.89)+
  facet_grid(~mod)
