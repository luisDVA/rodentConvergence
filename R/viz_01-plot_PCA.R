# plot PCA results 
# run PCA script first
library(dplyr)
library(readr)
library(tidyr)
library(tidytext)
library(FactoMineR)
library(ggplot2)
library(scico)
library(stringr)
library(patchwork)

# arrange by correlation
corrDatall <-
  dimscorrall %>%
  mutate(component=str_replace(component,"Dim.","PC")) %>% 
  mutate(trait = tidytext::reorder_within(trait, -correlation, within = component))

loadingsPltall <-
  ggplot(contrDatall, aes(x = trait, y = contribution, color = contribution)) +
  stat_summary(fun = "identity", position = "identity", geom = "crossbar", size = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 2) +
  scale_color_scico(palette = "imola", direction = -1) +
  ggthemes::theme_few(base_size = 12, base_family = "Laksaman") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()) +
  scale_x_reordered()+ggtitle("c)")+
  facet_grid(~component, scales = "free_x")
loadingsPltall

correlatPltall <-
  ggplot(corrDatall, aes(x = trait, y = correlation, color = correlation)) +
  stat_summary(fun = "identity", position = "identity", geom = "crossbar", size = 1, show.legend = FALSE) +
  scale_color_scico(palette = "tofino", direction = -1) +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 2) +
  ggthemes::theme_few(base_size = 12, base_family = "Laksaman") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank()) +
  labs(x = "measurement") +
  scale_x_reordered()+
  facet_grid(~component, scales = "free_x")
correlatPltall


traitsPCApltall <-
  pcssbdall %>% 
  ggplot() +
  geom_point(aes(Dim.1, Dim.2,group=sp),
             color = "black", pch = 21, size = 2) +
  #scale_fill_scico_d(palette = "nuuk") +
  geom_hline(yintercept = 0, size = 0.5, color = "dark grey", linetype = 3) +
  geom_vline(xintercept = 0, size = 0.5, color = "dark grey", linetype = 3) +
  coord_equal() +
  ggthemes::theme_few(base_size = 16, base_family = "Laksaman") 


traitsPCApltall

plotly::ggplotly(traitsPCApltall)
pcssbd

allpcsplot <- loadingsPltall+correlatPltall+plot_layout(nrow = 2)#+plot_annotation(title="c)")

allpcsplot
# all three
crplot
extplot
allpcsplot

(crplot|extplot)/allpcsplot
ggsave(filename = "2021_conv/figs/PCAall.png",width = 12,height = 7,units = "in",dpi = 300)
