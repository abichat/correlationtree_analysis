library(tidyverse)
library(cowplot)
library(ggtree)
library(scales)

tree_phy <- readRDS("real_datasets/chlamydiae/chlamydiae-tree_phy.rds")
tree_cor <- readRDS("real_datasets/chlamydiae/chlamydiae-tree_cor.rds")
tbl_pvalues <- readRDS("real_datasets/chlamydiae/chlamydiae-pvalues.rds")
df_abund_cor <- readRDS("real_datasets/chlamydiae/chlamydiae-abund_cor.rds")

mytheme <-
  theme_minimal() +
  theme(strip.text = element_text(size = 10),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) 

## Trees 

p_phy_phy <-
  tree_phy %>%
  ggtree(branch.length="none", color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(phy)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = phy < 0.1), vjust = -0.3, hjust = 1, 
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1) 

p_phy_cor <-
  tree_phy %>% 
  ggtree(branch.length="none", color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(cor)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = cor < 0.1), vjust = -0.3,
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1) +
  scale_x_reverse()

p_cor_phy <-
  tree_cor %>% 
  ggtree(color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(phy)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = phy < 0.1), vjust = -0.3, hjust = 1,
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1)

p_cor_cor <-
  tree_cor %>% 
  ggtree(color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(cor)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = cor < 0.1), vjust = -0.3,
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1) +
  scale_x_reverse()


p_trees <-
  plot_grid(p_phy_phy, p_phy_cor, NULL, NULL, p_cor_phy, p_cor_cor, 
            ncol = 2, rel_heights = c(1, 0.1, 1),
            labels = c("A", "B", "", "", "C", "D"), 
            label_x = c(0, 0.9, 0, 0, 0, 0.9), label_y = c(1, 1, 0, 0, 0.12, 0.12))

## Boxplots 

p_abund1 <- 
  df_abund_cor %>% 
  filter(OTU == 547579) %>% 
  ggplot() +
  aes(x = SampleType, y = Abundance) +
  geom_boxplot() +
  facet_wrap(~ OTU, scales = "free_y") +
  scale_y_sqrt(breaks = pretty_breaks(4)) +
  labs(x = NULL) +
  mytheme

p_abund2 <- 
  df_abund_cor %>% 
  filter(OTU != 547579) %>% 
  ggplot() +
  aes(x = SampleType, y = Abundance) +
  geom_boxplot() +
  facet_wrap(~ OTU, scales = "free_y") +
  scale_y_sqrt(breaks = pretty_breaks(4)) +
  labs(x = NULL, y = NULL) +
  mytheme

## Full plot 

p_trees <-
  plot_grid(p_phy_phy, p_phy_cor, NULL, NULL, p_cor_phy, p_cor_cor, 
            ncol = 2, rel_heights = c(1, 0.07, 1),
            labels = c("A", "B", "", "", "C", "D"), label_size = 14,
            label_x = c(0, 0.9, 0, 0, 0, 0.9), label_y = c(1, 1, 0, 0, 0.12, 0.12))

ggdraw(p_trees) +
  draw_plot(p_abund1, x = 0.01, y = 0.32, width = 0.3, height = 0.27) +
  draw_plot(p_abund2, x = 0.69, y = 0.32, width = 0.3, height = 0.27) +
  draw_plot_label(label = c("E", "F"), size = 14,
                  x = c(0.3, 0.66), y = c(0.55, 0.55)) +
  draw_plot_label("*", x = 0.59, y = c(0.442, 0.676), colour = "red")

ggsave("figures/figure_6.png", width = 15, height = 12, units = "cm")

