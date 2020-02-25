library(tidyverse)
library(cowplot)

df_4_boxplots <- readRDS("real_datasets/chaillou/chaillou-df_4_boxplots.rds")

bp_cor <-
  df_4_boxplots %>% 
  filter(Detected == "Correlation") %>% 
  ggplot() +
  aes(x = Env, y = Count) +
  geom_boxplot(aes(fill = Env)) +
  facet_wrap(vars(OTU), ncol = 6, scales = "free_y") +
  scale_fill_brewer(type = "qual") +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = NULL, fill = "Food type") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = alpha("blue", 0.5)), 
        axis.text.x = element_blank(),
        legend.position = "none")

legend_abund <- get_legend(bp_cor + theme(legend.position = "bottom"))

bp_phy <-
  df_4_boxplots %>% 
  filter(Detected == "Phylogeny") %>% 
  ggplot() +
  aes(x = Env, y = Count) +
  geom_boxplot(aes(fill = Env)) +
  facet_wrap(vars(OTU), ncol = 6, scales = "free_y") +
  scale_fill_brewer(type = "qual") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = alpha("firebrick", 0.5)),
        axis.text.x = element_blank(),
        legend.position = "none")

plot_grid(
  plot_grid(
    bp_cor, bp_phy, 
    ncol = 1, rel_heights = c(2, 1), 
    align = "v", axis = "l"),
  legend_abund, 
  ncol = 1, rel_heights = c(3, 0.1))

ggsave("figures/figure_s5.png", width = 15, height = 12, units = "cm")
