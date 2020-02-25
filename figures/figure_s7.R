library(tidyverse)
library(cowplot)
library(ggtree)

subtree_cor <- readRDS("real_datasets/chaillou/chaillou-tree_cob_sub.rds")
df_4_boxplots <- readRDS("real_datasets/chaillou/chaillou-df_4_boxplots.rds")
df_node_pv <- readRDS("real_datasets/chaillou/chaillou-df_node_pv.rds")

palette_detected <- setNames(c("blue", "purple", "firebrick", "black"), 
                             c("Correlation", "Both", "Phylogeny", "None"))

n_y <- 0.08
n_x <- 0.005
hj <- 0
col_lm <- "black"
col_kw <- "grey60"

p_subtree <-
  ggtree(subtree_cor) %<+% 
  df_node_pv +
  geom_point() +
  geom_tiplab(aes(color = detected), nudge_y = 0.1, hjust = 1.1) +
  geom_nodelab(aes(label = p_lm), color = col_lm, nudge_x = n_x, nudge_y = n_y, hjust = hj) +
  geom_nodelab(aes(label = p_kw), color = col_kw, nudge_x = n_x,  nudge_y = -n_y, hjust = hj) +
  geom_tiplab(aes(label = p_lm),  color = col_lm, offset = n_x,  nudge_y = n_y, hjust = hj) +
  geom_tiplab(aes(label = p_kw),  color = col_kw, offset = n_x,   nudge_y = -n_y, hjust =  hj) + 
  scale_color_manual(values = palette_detected) +
  xlim(0, 0.35)

bp_subtree <-
  df_4_boxplots %>% 
  filter(OTU %in% subtree_cor$tip.label) %>% 
  mutate(OTU = fct_relevel(OTU, "otu_00656", "otu_00519", "otu_01495", "otu_00516")) %>% 
  ggplot() +
  aes(x = Env, y = Count) +
  geom_boxplot(aes(fill = Env)) +
  facet_wrap(vars(OTU), ncol = 1, scales = "free_y") +
  scale_fill_brewer(type = "qual") +
  labs(x = NULL, y = NULL, fill = "Food type") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = alpha("grey80", 0.5)), 
        axis.text.x = element_blank(),
        legend.position = "right")

plot_grid(p_subtree, bp_subtree, ncol = 2)

ggsave("figures/figure_s7.png", width = 15, height = 15, units = "cm")
