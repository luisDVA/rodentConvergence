# size correction
library(readr)
library(dplyr)

# all species
meansmasses <- read_csv("2021_conv/data/meansMasses265_rdy.csv") 
# subset by measurment type
CRdat <- meansmasses %>% select(sp,CBL:ACP)
EXTdat <- meansmasses %>% select(sp,HB,T:UM)


# ACP separately
ACPcol <- CRdat %>% select(sp,ACP)
CRmorphlin <- CRdat %>% select(-ACP) 
ACPcol$ACP <- log(ACPcol$ACP)

# log shape ratio - craniodental & mandibular
crmorph_lsr <- log(
  CRmorphlin[,2:7]/ apply(CRmorphlin[,2:7], 1, prod)^(1/ncol(CRmorphlin[,2:7]))
  )
CRmorphdatf <- bind_cols(ACPcol,crmorph_lsr)


# log shape ratio - external
extmorph_lsr <- log(
  EXTdat[,2:8]/ apply(EXTdat[,2:8], 1, prod)^(1/ncol(EXTdat[,2:8]))
)

lsrall <- bind_cols(CRmorphdatf,extmorph_lsr)
lsrall <- lsrall %>% mutate(across(where(is.numeric),scale)) %>% 
  mutate(across(where(is.numeric),as.vector))
write_csv(lsrall,file = "2021_conv/data/logsize_rdy.csv")
