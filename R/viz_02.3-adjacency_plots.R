# Adjacency plots for network data capturing convergent regimes
library(purrr)
library(forcats)
library(paletteer)
library(patchwork)
library(ggstar)
library(RcppAlgos)

# Run 

plot_adj_tiles <- function(grphobject) {
  # order from graph
  adjorder <-
    grphobject %>%
    activate(nodes) %>%
    get.adjacency() %>%
    dimnames() %>%
    {
      .[[1]]
    } %>%
    {
      str_replace(., "chebezi", "ruschii")
    }

  datlongPCs <-
    grphobject %>%
    activate(edges) %>%
    mutate(
      to_name = .N()$name[to],
      from_name = .N()$name[from]
    ) %>%
    as_tibble() %>%
    select(from = from_name, to = to_name, ntrees, w, reg = reg.x) %>%
    mutate(w = round(w, 2)) %>%
    mutate(across(c(from, to), ~ str_replace(., "chebezi", "ruschii")))

  leglabs <-
    datlongPCs %>%
    group_by(reg) %>%
    distinct(reg, w) %>%
    summarize(reglab = paste0("Regime ", reg, ": ", "w =", w))

  datlongPCs <- left_join(datlongPCs, leglabs)


  datlongPCs <- datlongPCs %>% mutate(
    to = fct_drop(factor(to, levels = adjorder)),
    from = fct_drop(factor(from, levels = rev(adjorder)))
  )

  ghostdat <-
    tidyr::expand(datlongPCs, to, from) %>%
    select(to, from) %>%
    anti_join(., datlongPCs)


  theme_adjtiles <- function() {
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      text = element_text(family = "Laksaman"),
      panel.background = element_blank(),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(size = 16),
      axis.title = element_blank(),
      axis.text = element_text(face = "italic", color = "black", size = 9),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  }

  ggplot(data = datlongPCs) +
    geom_point(
      color = "black",
      aes(
        x = from, y = to, fill = ntrees,
        shape = reglab
      ), stroke = 0.3, size = 3.2
    ) +
    geom_point(
      data = ghostdat, color = "black",
      aes(x = from, y = to), pch = 45, size = 4
    ) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_paletteer_binned("scico::devon", direction = -1) +
    guides(
      shape = guide_legend(
        title.position = "bottom", nrow = 2, order = 1, title = "",
        override.aes = list(stroke = 1.4)
      ),
      fill = guide_colorbar(
        title.position = "bottom", title = "# of trees",
        frame.colour = "black",
        label.theme = element_text(size = 6),
        barheight = 0.5, barwidth = 4.2
      )
    ) +
    theme_adjtiles() +
    scale_x_discrete(expand = c(0.02, 0.02)) #+
  # scale_y_discrete(expand = c(0.02, 0.02))
}

plot_adj_tiles_wb <- function(grphobject) {
  # order from graph
  adjorder <-
    grphobject %>%
    activate(nodes) %>%
    get.adjacency() %>%
    dimnames() %>%
    {
      .[[1]]
    } %>%
    {
      str_replace(., "chebezi", "ruschii")
    }

  datlongPCs <-
    grphobject %>%
    activate(edges) %>%
    mutate(
      to_name = .N()$name[to],
      from_name = .N()$name[from]
    ) %>%
    as_tibble() %>%
    select(from = from_name, to = to_name, w, reg) %>%
    mutate(w = round(w, 2)) %>%
    mutate(across(c(from, to), ~ str_replace(., "chebezi", "ruschii")))

  leglabs <-
    datlongPCs %>%
    group_by(reg) %>%
    distinct(reg, w) %>%
    summarize(reglab = paste0("Regime ", reg, ": ", "w =", w))



  datlongPCs <- left_join(datlongPCs, leglabs)


  datlongPCs <- datlongPCs %>% mutate(
    to = fct_drop(factor(to, levels = adjorder)),
    from = fct_drop(factor(from, levels = rev(adjorder)))
  )

  ghostdat <-
    tidyr::expand(datlongPCs, to, from) %>%
    select(to, from) %>%
    anti_join(., datlongPCs)

  ghostnum <- ghostdat %>% mutate(across(everything(), as.numeric))

  selective_jitter <- function(x, # x = x co-ordinate
                               y, # y = y co-ordinate
                               g # g = group
  ) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    a <- cbind(x, y)
    a[duplicated(a)] <- jitter(a[duplicated(a)], amount = .12) # amount could be made a parameter

    final <- cbind(a, g)
    return(final)
  }

  newcoords <- as_tibble(selective_jitter(datlongPCs$from, datlongPCs$to, datlongPCs$reg))

  datlongPCs <-
    bind_cols(datlongPCs, newcoords)

  theme_adjtiles <- function() {
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_blank(),
      text = element_text(family = "Laksaman"),
      panel.background = element_blank(),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(size = 16),
      axis.title = element_blank(),
      axis.text = element_text(face = "italic", color = "black", size = 9),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  }

  frlabs <- datlongPCs$to %>%
    unique() %>%
    as.character()
  tolabs <- datlongPCs$from %>%
    unique() %>%
    as.character()

  ggplot(data = datlongPCs) +
    geom_star(
      color = "black",
      aes(
        x = x, y = y, fill = reglab,
        starshape = reglab
      ), size = 3.2
    ) +
    geom_point(
      data = ghostnum, color = "black",
      aes(x = from, y = to), pch = 45, size = 4
    ) +
    scale_starshape_manual(values = c(15, 13, 28, 11, 30, 5)) +
    scale_fill_scico_d(palette = "lapaz", direction = -1) +
    # guides(fill = guide_legend(title.position="bottom",nrow=2,order=1,title=""))+
    theme_adjtiles() +
    scale_y_continuous(breaks = 1:length(frlabs), labels = frlabs) +
    scale_x_continuous(
      breaks = 1:length(tolabs), labels = tolabs,
      expand = c(0.02, 0.02)
    )
}

plot_adj_tiles_flip <- function(grphobject) {
  # order from graph
  adjorder <-
    grphobject %>%
    activate(nodes) %>%
    get.adjacency() %>%
    dimnames() %>%
    {
      .[[1]]
    } %>%
    {
      str_replace(., "chebezi", "ruschii")
    }

  datlongPCs <-
    grphobject %>%
    activate(edges) %>%
    mutate(
      to_name = .N()$name[to],
      from_name = .N()$name[from]
    ) %>%
    as_tibble() %>%
    select(from = from_name, to = to_name, ntrees, w, reg = reg.x) %>%
    mutate(w = round(w, 2)) %>%
    mutate(across(c(from, to), ~ str_replace(., "chebezi", "ruschii")))

  leglabs <-
    datlongPCs %>%
    group_by(reg) %>%
    distinct(reg, w) %>%
    summarize(reglab = paste0("Regime ", reg, ": ", "w =", w))

  datlongPCs <- left_join(datlongPCs, leglabs)


  datlongPCs <- datlongPCs %>% mutate(
    to = fct_drop(factor(to, levels = adjorder)),
    from = fct_drop(factor(from, levels = rev(adjorder)))
  )

  ghostdat <-
    tidyr::expand(datlongPCs, to, from) %>%
    select(to, from) %>%
    anti_join(., datlongPCs)


  # datlongPC1cr %>% janitor::tabyl(reg,ntrees)
  # datlongPC1cr %>% janitor::tabyl(reg,w)

  theme_adjtiles <- function() {
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      text = element_text(family = "Laksaman"),
      panel.background = element_blank(),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(size = 16),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_blank(),
      axis.text = element_text(face = "italic", color = "black", size = 9),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  }

  ggplot(data = datlongPCs) +
    geom_point(
      color = "black",
      aes(
        x = from, y = to, fill = ntrees,
        shape = reglab
      ), stroke = 0.3, size = 3.2
    ) +
    geom_point(
      data = ghostdat, color = "black",
      aes(x = from, y = to), pch = 45, size = 4
    ) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_paletteer_binned("scico::devon", direction = -1) +
    guides(
      shape = guide_legend(
        title.position = "bottom", nrow = 2, title = "", order = 1,
        override.aes = list(stroke = 1.4)
      ),
      fill = guide_colorbar(
        title.position = "bottom", title = "# of trees",
        frame.colour = "black",
        label.theme = element_text(size = 6),
        barheight = 0.5, barwidth = 4.2
      )
    ) +
    theme_adjtiles() +
    scale_x_discrete(expand = c(0.02, 0.02)) +
    coord_flip()
  # scale_y_discrete(expand = c(0.02, 0.02))
}

plot_adj_tiles_small <- function(grphobject) {
  # order from graph
  adjorder <-
    grphobject %>%
    activate(nodes) %>%
    get.adjacency() %>%
    dimnames() %>%
    {
      .[[1]]
    }

  datlongPCs <-
    grphobject %>%
    activate(edges) %>%
    mutate(
      to_name = .N()$name[to],
      from_name = .N()$name[from]
    ) %>%
    as_tibble() %>%
    select(from = from_name, to = to_name, ntrees, w, reg = reg.x) %>%
    mutate(w = round(w, 2)) %>%
    mutate(across(c(from, to), ~ str_replace(., "chebezi", "ruschii")))

  leglabs <-
    datlongPCs %>%
    group_by(reg) %>%
    distinct(reg, w) %>%
    summarize(reglab = paste0("Regime ", reg, ": ", "w =", w))

  datlongPCs <- left_join(datlongPCs, leglabs)


  datlongPCs <- datlongPCs %>% mutate(
    to = fct_drop(factor(to, levels = adjorder)),
    from = fct_drop(factor(from, levels = rev(adjorder)))
  )

  ghostdat <-
    tidyr::expand(datlongPCs, to, from) %>%
    select(to, from) %>%
    anti_join(., datlongPCs)


  theme_adjtiles <- function() {
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      text = element_text(family = "At"),
      panel.background = element_blank(),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(size = 16),
      axis.title = element_blank(),
      axis.text = element_text(face = "italic", color = "black", size = 9),
      axis.text.x = element_text(hjust = 0.5)
    )
  }

  ggplot(data = datlongPCs) +
    geom_point(
      color = "black",
      aes(
        x = from, y = to, fill = ntrees,
        shape = reglab
      ), stroke = 0.3, size = 3.2
    ) +
    geom_point(
      data = ghostdat, color = "black",
      aes(x = from, y = to), pch = 45, size = 4
    ) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_paletteer_binned("scico::devon", direction = -1) +
    guides(
      shape = guide_legend(
        title.position = "bottom", nrow = 2, order = 1, title = "",
        override.aes = list(stroke = 1.4)
      ),
      fill = guide_colorbar(
        title.position = "bottom", title = "# of trees",
        frame.colour = "black",
        label.theme = element_text(size = 6),
        barheight = 0.5, barwidth = 4.2
      )
    ) +
    theme_adjtiles() +
    scale_x_discrete(expand = c(-0.005, 0)) #+
  # scale_y_discrete(expand = c(0, 0))
}


####
plot_adj_tiles(PC1regscrGrph) + PC1crregsPhylo +
  plot_layout(widths = c(0.5, 1.5))

# plot craniodental ----
plot_adj_tiles_flip(PC1regscrGrph) +
  plot_adj_tiles(PC2regscrGrph) +
  plot_adj_tiles(PC3regscrGrph) +
  plot_layout(widths = c(0.8, 1, 1)) + plot_annotation(tag_levels = "A")
ggsave("2021_conv/figs/new/PC1-3crregsAdj.png",
  device = agg_png,
  height = 7.3, width = 14.7, units = "in", dpi = 400
)

# plot external data ----
plot_adj_tiles_single <- function(grphobject) {
  # order from graph
  adjorder <-
    grphobject %>%
    activate(nodes) %>%
    get.adjacency() %>%
    dimnames() %>%
    {
      .[[1]]
    } %>%
    {
      str_replace(., "chebezi", "ruschii")
    }

  datlongPCs <-
    grphobject %>%
    activate(edges) %>%
    mutate(
      to_name = .N()$name[to],
      from_name = .N()$name[from]
    ) %>%
    as_tibble() %>%
    select(from = from_name, to = to_name, ntrees, w, reg = reg.x) %>%
    mutate(w = round(w, 2)) %>%
    mutate(across(c(from, to), ~ str_replace(., "chebezi", "ruschii")))

  leglabs <-
    datlongPCs %>%
    group_by(reg) %>%
    distinct(reg, w) %>%
    summarize(reglab = paste0("Regime ", reg, ": ", "w =", w))

  datlongPCs <- left_join(datlongPCs, leglabs)


  datlongPCs <- datlongPCs %>% mutate(
    to = fct_drop(factor(to, levels = adjorder)),
    from = fct_drop(factor(from, levels = rev(adjorder)))
  )

  ghostdat <-
    tidyr::expand(datlongPCs, to, from) %>%
    select(to, from) %>%
    anti_join(., datlongPCs)


  # datlongPC1cr %>% janitor::tabyl(reg,ntrees)
  # datlongPC1cr %>% janitor::tabyl(reg,w)

  theme_adjtiles <- function() {
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      text = element_text(family = "Laksaman"),
      panel.background = element_blank(),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(size = 16),
      axis.title = element_blank(),
      axis.text = element_text(face = "italic", color = "black", size = 9),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  }

  ggplot(data = datlongPCs) +
    geom_point(
      color = "black",
      aes(
        x = from, y = to, fill = ntrees,
        shape = reglab
      ), stroke = 0.3, size = 3.2
    ) +
    geom_point(
      data = ghostdat, color = "black",
      aes(x = from, y = to), pch = 45, size = 4
    ) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_scico(palette = "devon", direction = -1) +
    guides(
      shape = guide_legend(
        title.position = "bottom", nrow = 2, order = 1, title = "",
        override.aes = list(stroke = 1.4)
      ),
      fill = guide_colorbar(
        title.position = "right", title = "trees",
        frame.colour = "black",
        label.theme = element_text(size = 6),
        barheight = 0.5, barwidth = 1
      )
    ) +
    theme_adjtiles() +
    scale_x_discrete(expand = c(0.02, 0.02)) #+
  # scale_y_discrete(expand = c(0.02, 0.02))
}
plot_adj_tiles_3reg <- function(grphobject) {
  # order from graph
  adjorder <-
    grphobject %>%
    activate(nodes) %>%
    get.adjacency() %>%
    dimnames() %>%
    {
      .[[1]]
    } %>%
    {
      str_replace(., "chebezi", "ruschii")
    }

  datlongPCs <-
    grphobject %>%
    activate(edges) %>%
    mutate(
      to_name = .N()$name[to],
      from_name = .N()$name[from]
    ) %>%
    as_tibble() %>%
    select(from = from_name, to = to_name, ntrees, w, reg = reg.x) %>%
    mutate(w = round(w, 2)) %>%
    mutate(across(c(from, to), ~ str_replace(., "chebezi", "ruschii")))

  leglabs <-
    datlongPCs %>%
    group_by(reg) %>%
    distinct(reg, w) %>%
    summarize(reglab = paste0("Regime ", reg, ": ", "w =", w))

  datlongPCs <- left_join(datlongPCs, leglabs)


  datlongPCs <- datlongPCs %>% mutate(
    to = fct_drop(factor(to, levels = adjorder)),
    from = fct_drop(factor(from, levels = rev(adjorder)))
  )

  ghostdat <-
    tidyr::expand(datlongPCs, to, from) %>%
    select(to, from) %>%
    anti_join(., datlongPCs)


  # datlongPC1cr %>% janitor::tabyl(reg,ntrees)
  # datlongPC1cr %>% janitor::tabyl(reg,w)

  theme_adjtiles <- function() {
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      text = element_text(family = "Laksaman"),
      panel.background = element_blank(),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(size = 16),
      axis.title = element_blank(),
      axis.text = element_text(face = "italic", color = "black", size = 9),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  }

  ggplot(data = datlongPCs) +
    geom_point(
      color = "black",
      aes(
        x = from, y = to, fill = ntrees,
        shape = reglab
      ), stroke = 0.3, size = 3.2
    ) +
    geom_point(
      data = ghostdat, color = "black",
      aes(x = from, y = to), pch = 45, size = 4
    ) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_paletteer_binned("scico::devon", direction = -1) +
    guides(
      shape = guide_legend(
        title.position = "bottom", nrow = 3, order = 1, title = "",
        override.aes = list(stroke = 1.4)
      ),
      fill = guide_colorbar(
        title.position = "bottom", title = "# of trees",
        frame.colour = "black",
        label.theme = element_text(size = 6),
        barheight = 0.5, barwidth = 4.2
      )
    ) +
    theme_adjtiles() +
    scale_x_discrete(expand = c(0.02, 0.02)) #+
  # scale_y_discrete(expand = c(0.02, 0.02))
}

plot_adj_tiles_3reg(PC1regsextGrph) +
  plot_adj_tiles_single(PC2regsextGrph) +
  plot_adj_tiles_single(PC3regsextGrph) +
  # plot_layout(widths=c(1.1,0.9,0.9))+
  plot_annotation(tag_levels = "A")
ggsave("2021_conv/figs/new/PC1-3extregsAdj.png",
  device = agg_png,
  height = 7.3, width = 14.7, units = "in", dpi = 400
)
plot_adj_tiles_small(MultregscrGrph) + plot_layout(nrow = 2)

# plot multivariate regs ----
((plot_adj_tiles_small(MultregscrGrph) /
  plot_adj_tiles(MultregsextGrph)) |
  plot_adj_tiles_wb(wbregsGrph)) +
  plot_layout(widths = c(0.45, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("2021_conv/figs/new/MultregsAdj.png",
  device = agg_png,
  height = 7.3, width = 12.7, units = "in", dpi = 400
)
