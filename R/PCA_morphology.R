# PCA on morphology data
library(dplyr)
library(readr)
library(tidyr)
library(tidytext)
library(FactoMineR)
library(stringr)


# morphology data (size corrected)
allmorph <- read_csv("2021_conv/data/logsize_rdy.csv")

# subset craniondental and mandibular
crdat <- allmorph %>% select(ACP:HMC)
extdat <- allmorph %>% select(HB:UM)



# cranial -----------------------------------------------------------------
# calculate PCS and interpret
crdat_pca <- PCA(as.data.frame(crdat))
summary(crdat_pca)
# get descriptions for each dimension
dimsload <- dimdesc(crdat_pca)

# get correlations
dim1corr <- as_tibble(dimsload$Dim.1[[1]], rownames = "trait") %>% mutate(component = "Dim.1")
dim2corr <- as_tibble(dimsload$Dim.2[[1]], rownames = "trait") %>% mutate(component = "Dim.2")
dim3corr <- as_tibble(dimsload$Dim.3[[1]], rownames = "trait") %>% mutate(component = "Dim.3")
dimscorr <- bind_rows(dim1corr, dim2corr,dim3corr)
# variable contributions
contribs <- crdat_pca$var$contrib
# PC coordinates
pccooords <- crdat_pca$ind$coord

# together
pcssbd <- bind_cols(allmorph, as_tibble(pccooords))
pcssbd
# for export
pcssbd %>% dplyr::select(sp,contains("Dim")) -> pcsCr

# melt
contribslong <- contribs %>%
  as_tibble(rownames = "trait") %>%
  pivot_longer(
    cols = Dim.1:Dim.5, names_to = "component",
    values_to = "contribution")

# arrange by contribution
contrDat <-
  contribslong %>%
  filter(component == "Dim.1" | component == "Dim.2"|component == "Dim.3") %>%
  mutate(component=str_replace(component,"Dim.","PC")) %>% 
  mutate(trait = tidytext::reorder_within(trait, -contribution, within = component))

# arrange by correlation
corrDat <-
  dimscorr %>%
  mutate(trait = tidytext::reorder_within(trait, -correlation, within = component)) %>% 
  mutate(component=str_replace(component,"Dim.","PC")) 

loadingsPlt <-
  ggplot(contrDat, aes(x = trait, y = contribution, color = contribution)) +
  stat_summary(fun = "identity", position = "identity", geom = "crossbar", size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 2) +
  scale_color_scico(palette = "imola", direction = -1) +
  ggthemes::theme_few(base_size = 12, base_family = "Laksaman") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  scale_x_reordered()+ggtitle("a)")+
  facet_grid(~component, scales = "free_x")
loadingsPlt

correlatPlt <-
  ggplot(corrDat, aes(x = trait, y = correlation, color = correlation)) +
  stat_summary(fun = "identity", position = "identity", geom = "crossbar", size = 1, show.legend = FALSE) +
  scale_color_scico(palette = "tofino", direction = -1) +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 2) +
  ggthemes::theme_few(base_size = 12, base_family = "Laksaman") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank()) +
  scale_x_reordered()+
  labs(x = "measurement") +
  facet_grid(~component, scales = "free_x")
correlatPlt


traitsPCAplt <-
  pcssbd %>% 
  ggplot() +
  geom_point(aes(Dim.1, Dim.2,group=sp),
             color = "black", pch = 21, size = 2) +
  #scale_fill_scico_d(palette = "nuuk") +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 3) +
  geom_vline(xintercept = 0, size = 0.5, color = "dark grey", linetype = 3) +
  coord_equal() +
  ggthemes::theme_few(base_size = 16, base_family = "Laksaman") 
  
traitsPCAplt
pcssbd %>% 
  ggplot() +
  geom_point(aes(Dim.1, CBL,group=sp),
             color = "black", pch = 21, size = 1)+
geom_point(aes(Dim.1, LR,group=sp),
           color = "red", pch = 21, size = 1) +
  geom_point(aes(Dim.1, ACP,group=sp),
             color = "green", pch = 21, size = 2) 

plotly::ggplotly(traitsPCAplt)
pcssbd


# PC1
plot(pcssbd$ACP,pcssbd$Dim.1)


crplot <- loadingsPlt+correlatPlt+plot_layout(nrow=2)#+plot_annotation(title="a)")
crplot

# external ----------------------------------------------------------------
# calculate PCS and interpret
extdat_pca <- PCA(as.data.frame(extdat))
summary(extdat_pca)
# get descriptions for each dimension
dimsloadext <- dimdesc(extdat_pca)

# get correlations
dim1corrext <- as_tibble(dimsloadext$Dim.1[[1]], rownames = "trait") %>% mutate(component = "Dim.1")
dim2corrext <- as_tibble(dimsloadext$Dim.2[[1]], rownames = "trait") %>% mutate(component = "Dim.2")
dim3corrext <- as_tibble(dimsloadext$Dim.3[[1]], rownames = "trait") %>% mutate(component = "Dim.3")
dimscorrext <- bind_rows(dim1corrext, dim2corrext,dim3corrext)
# variable contributions
contribsext <- extdat_pca$var$contrib
# PC coordinates
pccooordsext <- extdat_pca$ind$coord

# together
pcssbdext <- bind_cols(allmorph, as_tibble(pccooordsext))

# for export
pcssbdext %>% dplyr::select(sp,contains("Dim")) -> pcsExt

# melt
contribslongext <- contribsext %>%
  as_tibble(rownames = "trait") %>%
  pivot_longer(
    cols = Dim.1:Dim.5, names_to = "component",
    values_to = "contribution")

# arrange by contribution
contrDatext <-
  contribslongext %>%
  filter(component == "Dim.1" | component == "Dim.2"|component == "Dim.3") %>%
  mutate(component=str_replace(component,"Dim.","PC")) %>% 
  mutate(trait = tidytext::reorder_within(trait, -contribution, within = component))

# arrange by correlation
corrDatext <-
  dimscorrext %>%
  mutate(component=str_replace(component,"Dim.","PC")) %>% 
  mutate(trait = tidytext::reorder_within(trait, -correlation, within = component))

loadingsPltext <-
  ggplot(contrDatext, aes(x = trait, y = contribution, color = contribution)) +
  stat_summary(fun = "identity", position = "identity", geom = "crossbar", size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 2) +
  scale_color_scico(palette = "imola", direction = -1) +
  ggthemes::theme_few(base_size = 12, base_family = "Laksaman") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  scale_x_reordered()+ggtitle("b)")+
  facet_grid(~component, scales = "free_x")
loadingsPltext

correlatPltext <-
  ggplot(corrDatext, aes(x = trait, y = correlation, color = correlation)) +
  stat_summary(fun = "identity", position = "identity", geom = "crossbar", size = 1, show.legend = FALSE) +
  scale_color_scico(palette = "tofino", direction = -1) +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 2) +
  ggthemes::theme_few(base_size = 12, base_family = "Laksaman") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank()) +
  labs(x = "measurement") +
  facet_grid(~component, scales = "free_x")+
  scale_x_reordered()
correlatPltext


traitsPCApltext <-
  pcssbdext %>% 
  ggplot() +
  geom_point(aes(Dim.2, Dim.3,group=sp),
             color = "black", pch = 21, size = 2) +
  #scale_fill_scico_d(palette = "nuuk") +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 3) +
  geom_vline(xintercept = 0, size = 0.5, color = "dark grey", linetype = 3) +
  coord_equal() +
  ggthemes::theme_few(base_size = 16, base_family = "Laksaman") 

traitsPCApltext

plotly::ggplotly(traitsPCApltext)
pcssbd


extplot <- loadingsPltext+correlatPltext+plot_layout(nrow=2)#+plot_annotation(title="b)")
extplot

# all traits ---- 

# calculate PCS and interpret
alldat_pca <- PCA(as.data.frame(allmorph[,2:15]))
summary(alldat_pca)
# get descriptions for each dimension
dimsloadall <- dimdesc(alldat_pca,1:5)
dimsloadall
# get correlations
dim1corrall <- as_tibble(dimsloadall$Dim.1[[1]], rownames = "trait") %>% mutate(component = "Dim.1")
dim2corrall <- as_tibble(dimsloadall$Dim.2[[1]], rownames = "trait") %>% mutate(component = "Dim.2")
dim3corrall <- as_tibble(dimsloadall$Dim.3[[1]], rownames = "trait") %>% mutate(component = "Dim.3")
dim4corrall <- as_tibble(dimsloadall$Dim.4[[1]], rownames = "trait") %>% mutate(component = "Dim.4")
dim5corrall <- as_tibble(dimsloadall$Dim.5[[1]], rownames = "trait") %>% mutate(component = "Dim.5")

dimscorrall <- bind_rows(dim1corrall, dim2corrall,dim3corrall,dim4corrall,dim5corrall)
# variable contributions
contribsall <- alldat_pca$var$contrib
# PC coordinates
pccooordsall <- alldat_pca$ind$coord
pccooordsall
# together
pcssbdall <- bind_cols(allmorph, as_tibble(pccooordsall))
pcssbdall

# for export
pcssbdall %>% dplyr::select(sp,contains("Dim")) -> pcssbdall

# melt
contribslongall <- contribsall %>%
  as_tibble(rownames = "trait") %>%
  pivot_longer(
    cols = Dim.1:Dim.5, names_to = "component",
    values_to = "contribution")

# arrange by contribution
contrDatall <-
  contribslongall %>%
  filter(component == "Dim.1" | component == "Dim.2"|component == "Dim.3"|component=="Dim.4"|component=="Dim.5") %>%
  mutate(component=str_replace(component,"Dim.","PC")) %>% 
  mutate(trait = tidytext::reorder_within(trait, -contribution, within = component))


# export PCS
write_csv(pcsCr,"2021_conv/data/crPCs_fm.csv")
write_csv(pcsExt,"2021_conv/data/extPCs_fm.csv")
write_csv(pcssbdall,"2021_conv/data/allPCs_fm.csv")

