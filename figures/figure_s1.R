library(tidyverse)
library(cowplot)

source("figures/theme.R")
mytheme <- theme_minimal() + 
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 4.5),
        legend.position = "none")

tree_labels_chlamydia <- readRDS("forests/chlamydiae/chlamydiae-tree-labels.rds")
dist_bhv_chlamydia <- readRDS("forests/chlamydiae/chlamydiae-dist-bhv.rds")
pcoa_bhv_chlamydia <- readRDS("forests/chlamydiae/chlamydiae-pcoa-bhv.rds")

dist_df_bhv_chlamydia <- tibble(Distance = dist_bhv_chlamydia[1:(length(tree_labels_chlamydia) - 1)], 
                            Type = tree_labels_chlamydia[-1])
dist_tax_bhv_chlamydia <- filter(dist_df_bhv_chlamydia, Type == "Phylogeny")$Distance

p_bp_chlamydiae <-
  dist_df_bhv_chlamydia %>% 
  filter(!Type %in% c("Correlation", "Phylogeny")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_tax_bhv_chlamydia, color = color_values["Phylogeny"]) +
  geom_text(x = 1, y = dist_tax_bhv_chlamydia, vjust = 1.2, size = 2,
            color = color_values["Phylogeny"], label = "Phylogeny") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  scale_x_discrete(labels = c("Bootstrap", "Random\nCorrelation", "Random\nPhylogeny")) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme +
  theme(axis.text.x = element_text(size = 6))

p_pcoa_chlamydiae <-
  pcoa_bhv_chlamydia$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_labels_chlamydia) %>% 
  ggplot() +
  aes(Axis.1, Axis.2, 
      color = Type, shape = Type, size = Type, alpha = Type) +
  geom_point() +
  scale_color_manual(values = color_values) +
  scale_alpha_manual(values = alpha_values) +
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("Axis 1 (",round(pcoa_bhv_chlamydia$values$Relative_eig[1]*100, 2), " %)"),
       y = paste0("Axis 2 (",round(pcoa_bhv_chlamydia$values$Relative_eig[2]*100, 2), " %)")) +
  mytheme

plot_grid(p_pcoa_chlamydiae, p_bp_chlamydiae)

ggsave("figures/figure_s1.png", width = 15, height = 8, units = "cm")
