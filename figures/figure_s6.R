library(tidyverse)
library(cowplot)
library(ggtree)

tree_cor <- readRDS("real_datasets/chaillou/chaillou-tree_cor.rds")
tree_phy <- readRDS("real_datasets/chaillou/chaillou-tree_phy.rds")
df_pvalues_filtered <- readRDS("real_datasets/chaillou/chaillou-df_pvalues_filtered.rds")

palette_detected <- setNames(c("blue", "purple", "firebrick", "black"), 
                             c("Correlation", "Both", "Phylogeny", "None"))

p1 <-
  tree_cor %>%
  ggtree(branch.length = "none", color = "grey30") %<+%
  df_pvalues_filtered +
  geom_hilight(node = 163, fill = "green", alpha = 0.5) +
  geom_tippoint(aes(size = -log10(cor), color = Detected), alpha = 0.5) +
  geom_tiplab(aes(label = OTU_displayed, color = Detected), vjust = 0, hjust = 1, 
              size = 2.5, fontface = "bold", show.legend = FALSE) +
  guides(size = FALSE) +
  scale_color_manual(values = palette_detected, name = "Detected by")

legend_tree <- get_legend(p1 + theme(legend.position = "bottom"))

p2 <-
  tree_phy %>% 
  ggtree(color = "grey30") %<+%
  df_pvalues_filtered +
  geom_tippoint(aes(size = -log10(phy), color = Detected), alpha = 0.5) +
  geom_tiplab(aes(label = OTU_displayed, color = Detected), vjust = 0, hjust = 0,
              size = 2.5, fontface = "bold") +
  scale_color_manual(values = palette_detected) + 
  theme(legend.position = "none") + 
  scale_x_reverse()

plot_grid(
  plot_grid(p1, p2, 
            ncol = 2, 
            labels = c("Correlation tree", "Phylogeny"), 
            label_x = c(0, 0.3), 
            label_y = c(1, 1)),
  legend_tree, 
  ncol = 1, rel_heights = c(2, 0.06))

ggsave("figures/figure_s6.png", width = 15, height = 22, units = "cm")
