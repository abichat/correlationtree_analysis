library(tidyverse)
library(cowplot)

source("figures/theme.R")
mytheme <- theme_minimal() + 
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 4.5),
        legend.position = "none")

# Ravel

tree_labels_ravel <- readRDS("forests/ravel/ravel-tree-labels.rds")
dist_rf_ravel <- readRDS("forests/ravel/ravel-dist-rf.rds")
pcoa_rf_ravel <- readRDS("forests/ravel/ravel-pcoa-rf.rds")

dist_df_rf_ravel <- tibble(Distance = dist_rf_ravel[1:(length(tree_labels_ravel) - 1)], 
                            Type = tree_labels_ravel[-1])
dist_tax_rf_ravel <- filter(dist_df_rf_ravel, Type == "Taxonomy")$Distance

p_bp_ravel <-
  dist_df_rf_ravel %>% 
  filter(!Type %in% c("Correlation", "Taxonomy")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_tax_rf_ravel, color = color_values["Taxonomy"]) +
  geom_text(x = 1, y = dist_tax_rf_ravel, vjust = 1.2, size = 2,
            color = color_values["Taxonomy"], label = "Taxonomy") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  scale_x_discrete(labels = c("Bootstrap", "Random\nCorrelation", "Random\nTaxonomy")) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme +
  theme(axis.text.x = element_text(size = 6))

p_pcoa_ravel <-
  pcoa_rf_ravel$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_labels_ravel) %>% 
  ggplot() +
  aes(Axis.1, Axis.2, 
      color = Type, shape = Type, size = Type, alpha = Type) +
  geom_point() +
  scale_color_manual(values = color_values) +
  scale_alpha_manual(values = alpha_values) +
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("Axis 1 (",round(pcoa_rf_ravel$values$Relative_eig[1]*100, 2), " %)"),
       y = paste0("Axis 2 (",round(pcoa_rf_ravel$values$Relative_eig[2]*100, 2), " %)")) +
  mytheme

# Zeller

tree_labels_zeller <- readRDS("forests/zeller/zeller-tree-labels.rds")
dist_rf_zeller <- readRDS("forests/zeller/zeller-dist-rf.rds")
pcoa_rf_zeller <- readRDS("forests/zeller/zeller-pcoa-rf.rds")

dist_df_rf_zeller <- tibble(Distance = dist_rf_zeller[1:(length(tree_labels_zeller) - 1)], 
                             Type = tree_labels_zeller[-1])
dist_tax_rf_zeller <- filter(dist_df_rf_zeller, Type == "Taxonomy")$Distance

p_bp_zeller <-
  dist_df_rf_zeller %>% 
  filter(!Type %in% c("Correlation", "Taxonomy")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_tax_rf_zeller, color = color_values["Taxonomy"]) +
  geom_text(x = 1, y = dist_tax_rf_zeller, vjust = 1.2, size = 2,
            color = color_values["Taxonomy"], label = "Taxonomy") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  scale_x_discrete(labels = c("Bootstrap", "Random\nCorrelation", "Random\nTaxonomy")) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme +
  theme(axis.text.x = element_text(size = 6))

p_pcoa_zeller <-
  pcoa_rf_zeller$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_labels_zeller) %>% 
  ggplot() +
  aes(Axis.1, Axis.2, 
      color = Type, shape = Type, size = Type, alpha = Type) +
  geom_point() +
  scale_color_manual(values = color_values) +
  scale_alpha_manual(values = alpha_values) +
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("Axis 1 (",round(pcoa_rf_zeller$values$Relative_eig[1]*100, 2), " %)"),
       y = paste0("Axis 2 (",round(pcoa_rf_zeller$values$Relative_eig[2]*100, 2), " %)")) +
  mytheme

# Zeller

tree_labels_chaillou <- readRDS("forests/chaillou/chaillou-tree-labels.rds")
dist_rf_chaillou <- readRDS("forests/chaillou/chaillou-dist-rf.rds")
pcoa_rf_chaillou <- readRDS("forests/chaillou/chaillou-pcoa-rf.rds")

dist_df_rf_chaillou <- tibble(Distance = dist_rf_chaillou[1:(length(tree_labels_chaillou) - 1)], 
                               Type = tree_labels_chaillou[-1])
dist_tax_rf_chaillou <- filter(dist_df_rf_chaillou, Type == "Phylogeny")$Distance

p_bp_chaillou <-
  dist_df_rf_chaillou %>% 
  filter(!Type %in% c("Correlation", "Phylogeny")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_tax_rf_chaillou, color = color_values["Phylogeny"]) +
  geom_text(x = 1, y = dist_tax_rf_chaillou, vjust = 1.2, size = 2,
            color = color_values["Phylogeny"], label = "Phylogeny") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  scale_x_discrete(labels = c("Bootstrap", "Random\nCorrelation", "Random\nPhylogeny")) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme +
  theme(axis.text.x = element_text(size = 6))

p_pcoa_chaillou <-
  pcoa_rf_chaillou$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_labels_chaillou) %>% 
  ggplot() +
  aes(Axis.1, Axis.2, 
      color = Type, shape = Type, size = Type, alpha = Type) +
  geom_point() +
  scale_color_manual(values = color_values) +
  scale_alpha_manual(values = alpha_values) +
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("Axis 1 (",round(pcoa_rf_chaillou$values$Relative_eig[1]*100, 2), " %)"),
       y = paste0("Axis 2 (",round(pcoa_rf_chaillou$values$Relative_eig[2]*100, 2), " %)")) +
  mytheme

# All

plot_grid(
  NULL, NULL, NULL,
  p_bp_ravel, p_bp_zeller, p_bp_chaillou,
  p_pcoa_ravel, p_pcoa_zeller, p_pcoa_chaillou,
  ncol = 3,
  rel_heights = c(.1, 1, 1),
  rel_widths = c(1, 1, 1),
  labels = c("Ravel (Vagina)", "Zeller (Gut)", "Chaillou (8 Food Types)", 
             "", "", "", 
             "", "", ""), 
  label_x = c(0.2, 0.25, -0.05), label_y = 0.9, label_size = 8)

ggsave("figures/figure_s2.png", width = 15, height = 8, units = "cm")
