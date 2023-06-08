# SA data, Nov 2021
library(readr)
library(dplyr)

extdata <- read_csv("2021_conv/data/external_out.csv",locale = locale(encoding = "UTF-16"))
crdata <- read_csv("2021_conv/data/cranial_out.csv",locale = locale(encoding = "UTF-16"))

# for taxonomy
cwalk <- read_csv("2021_conv/data/specimens_crosswalk.csv") %>% select(-ID)

extjoined <- left_join(extdata,cwalk) %>% 
  select(species,new_sp,institution,specimenID,cat_no,sex,everything())

crjoined <- left_join(crdata,cwalk) %>% 
    select(species,new_sp,institution,specimenID,cat_no,sex,everything())

externalF <- extjoined %>% mutate(HBnew=LT-T,LTnew=HB+T) %>% 
  mutate(HB2=coalesce(HB,HBnew),LT2=coalesce(LT,LTnew)) %>%
  relocate(LT2,.before = LT) %>% relocate(HB2,.before = HB) %>% 
  select(species:sex,LT=LT2,HB=HB2,T:label_species)

extmeans <- externalF %>% group_by(new_sp) %>% summarise(across(LT:W,mean,na.rm=TRUE))
extmeans
crmeans <- crjoined %>% group_by(new_sp) %>% 
  summarise(across(CBL:ACP,mean,na.rm=TRUE)) %>% select(-LmD)
crmeans

meansmassesSA <- left_join(extmeans,crmeans)

# to disk
write_csv(meansmassesSA,"2021_conv/data/meansMassesSAnov2021.csv")

# read again
meansmassesSA <- read_csv("2021_conv/data/meansMassesSAnov2021.csv")
meansmasses17 <- read_csv("2021_conv/data/meansMasses2017.csv")

names(meansmassesSA)
names(meansmasses17)

meansmassesSA <- meansmassesSA %>% 
  mutate(dataProv="LD") %>% 
  rename(sp=new_sp,UM=Um,massGrams=W)

rodsnew <- bind_rows(meansmasses17,meansmassesSA)

rodsnew <- rodsnew %>% mutate(across(where(is.numeric),round,2))
# export
write_csv(rodsnew,"2021_conv/data/meansMasses265.csv")
