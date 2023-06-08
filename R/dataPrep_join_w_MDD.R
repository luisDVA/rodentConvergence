library(readr)
library(dplyr)

# specimens and taxonomy
mmasses <- read_csv("2021_conv/data/meansMasses265_rdy.csv") %>% 
  select(sp)
mdd <- read_csv("2021_conv/data/MDD/MDD_v1.9.1_6596species.csv") %>%
            filter(order=="RODENTIA")  
mdd <- mdd %>% select(sciName,order:genus)

taxhigher <- left_join(mmasses,mdd,by=c("sp"="sciName"))
write_csv(taxhigher,"joinedout.csv")
taxhigher %>% distinct(new_sp,.keep_all = TRUE)

spMDDnames <- read_csv("joinedout.csv")

spMDDnames <- spMDDnames %>% mutate(mdd_name=ifelse(is.na(mdd_name),sp,mdd_name))
spMDDnames
spMDDnamesTax <- left_join(spMDDnames,mdd,by=c("mdd_name"="sciName"))

spMDDnamesTax                      
write_csv(spMDDnamesTax,file = "2021_conv/data/spMDD2022taxonomy.csv")

