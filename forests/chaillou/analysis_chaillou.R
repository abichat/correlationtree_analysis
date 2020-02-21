library(correlationtree)
library(biomformat)
library(tidyverse)
library(phyloseq)
library(distory)
library(broom)
library(furrr)
library(ape)

# Parallelize the code
options("future.fork.enable" = TRUE)
plan(multiprocess)


#### Data ####

biom <- read_biom("forests/chaillou/chaillou.biom")

## Abundance

prevalence_min <- 0.05

df_abund <- 
  biom %>% 
  biom_data() %>% 
  as("matrix") %>% 
  as_tibble(rownames = "OTU") %>%
  gather(key = "sample", value = "abundance", -OTU) %>% 
  group_by(OTU) %>% 
  mutate(P = mean(abundance > 0)) %>% 
  filter(P > prevalence_min) %>% 
  select(-P) %>% 
  ungroup() %>% 
  spread(key = sample, value = abundance)


## Trees 

tree_cor <- correlation_tree(df_abund, method = "spearman")
mean_lineage <- mean_lineage_length(tree_cor)

tree_phy <- 
  read.tree("forests/chaillou/phytree_chaillou.nwk") %>% 
  prune_taxa(df_abund$OTU, .)

tree_phy$edge.length <- tree_phy$edge.length * mean_lineage / mean_lineage_length(tree_phy)


#### Forest ####

set.seed(42)

N_boot <- 100 
N_rand <- 100

## Bootstraps

correlation_tree <- possibly(correlation_tree, NULL)

trees_boot <- 
  N_boot %>% 
  rerun(sample_boot(df_abund)) %>% 
  future_map(correlation_tree, method = "spearman") %>%
  discard(is.null) %>% 
  reduce(c)

N_boot <- length(trees_boot)

## Random trees 

trees_rand_cor <- 
  N_rand %>% 
  rerun(shuffle_tiplabels(tree_cor)) %>% 
  reduce(c)

trees_rand_phy <-
  N_rand %>%
  rerun(shuffle_tiplabels(tree_phy)) %>%
  map(multi2di) %>%
  reduce(c)

## Aggregation

forest <- c(tree_cor, tree_phy, trees_boot, trees_rand_cor, trees_rand_phy)

tree_labels <- 
  factor(c("Correlation", "Phylogeny", 
           rep("Bootstrap", N_boot), 
           rep("Random Correlation", N_rand), 
           rep("Random Phylogeny", N_rand)),
         levels = c("Correlation", "Bootstrap", 
                    "Random Correlation", "Random Phylogeny", "Phylogeny"))

# saveRDS(tree_labels, "forests/chaillou/chaillou-tree-labels.rds")
# tree_labels <- readRDS("forests/chaillou/chaillou-tree-labels.rds")

#### Distances and models ####

# Pairwise distances, could take time
# dist_bhv <- future_dist_BHV(forest) # Billera-Holmes-Vogtmann
dist_bhv <- dist.multiPhylo(forest) # Billera-Holmes-Vogtmann
dist_rf <- dist.topo(unroot(forest)) # Robinson-Foulds

# saveRDS(dist_bhv, "forests/chaillou/chaillou-dist-bhv.rds")
# saveRDS(dist_rf,  "forests/chaillou/chaillou-dist-rf.rds")
# dist_bhv <- readRDS("forests/chaillou/chaillou-dist-bhv.rds")
# dist_rf <-  readRDS("forests/chaillou/chaillou-dist-rf.rds")

# PCoA
pcoa_bhv <- pcoa(dist_bhv) 
pcoa_rf <- pcoa(dist_rf) 

# saveRDS(pcoa_bhv, "forests/chaillou/chaillou-pcoa-bhv.rds")
# saveRDS(pcoa_rf,  "forests/chaillou/chaillou-pcoa-rf.rds")
# pcoa_bhv <- readRDS("forests/chaillou/chaillou-pcoa-bhv.rds")
# pcoa_rf <-  readRDS("forests/chaillou/chaillou-pcoa-rf.rds")

# Distances to correlation tree
dist_df_bhv <- tibble(Distance = dist_bhv[1:(length(tree_labels) - 1)], Type = tree_labels[-1])
dist_df_rf  <- tibble(Distance = dist_rf[1:(length(tree_labels) - 1)],  Type = tree_labels[-1])

# Linear models
lm_bhv <- lm(Distance ~ Type, data = dist_df_bhv)
lm_rf  <- lm(Distance ~ Type, data = dist_df_rf)

aov_bhv <- aov(Distance ~ Type, data = dist_df_bhv)
aov_rf  <- aov(Distance ~ Type, data = dist_df_rf)


#### Results ####

glance(lm_bhv)
glance(lm_rf)

as_tibble(TukeyHSD(aov_bhv, "Type")$Type, rownames = "X")
as_tibble(TukeyHSD(aov_rf, "Type")$Type, rownames = "X")


#### Plots ####

## Themes 

source("figures/theme.R")

## Boxplots

dist_phy_bhv <- filter(dist_df_bhv, Type == "Phylogeny")$Distance
dist_phy_rf  <- filter(dist_df_rf, Type == "Phylogeny")$Distance

dist_df_bhv %>% 
  filter(!Type %in% c("Correlation", "Phylogeny")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_phy_bhv, color = color_values["Phylogeny"]) +
  geom_text(x = 1, y = dist_phy_bhv, vjust = 1.2, size = 4,
            color = color_values["Phylogeny"], label = "Phylogeny") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme

ggsave("forests/chaillou/chaillou-boxplot-bhv.png", width = 7.5, height = 5, dpi = "retina")

dist_df_rf %>% 
  filter(!Type %in% c("Correlation", "Phylogeny")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_phy_rf, color = color_values["Phylogeny"]) +
  geom_text(x = 1, y = dist_phy_rf, vjust = 1.2, size = 4,
            color = color_values["Phylogeny"], label = "Phylogeny") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme

ggsave("forests/chaillou/chaillou-boxplot-rf.png", width = 7.5, height = 5, dpi = "retina")


## PCoAs

pcoa_bhv$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_labels) %>% 
  ggplot() +
  aes(Axis.1, Axis.2, 
      color = Type, shape = Type, size = Type, alpha = Type) +
  geom_point() +
  scale_color_manual(values = color_values) +
  scale_alpha_manual(values = alpha_values) +
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("Axis 1 (",round(pcoa_bhv$values$Relative_eig[1]*100, 2), " %)"),
       y = paste0("Axis 2 (",round(pcoa_bhv$values$Relative_eig[2]*100, 2), " %)")) +
  mytheme

ggsave("forests/chaillou/chaillou-pcoa-bhv.png", width = 7.5, height = 5, dpi = "retina")

pcoa_rf$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_labels) %>% 
  ggplot() +
  aes(Axis.1, Axis.2, 
      color = Type, shape = Type, size = Type, alpha = Type) +
  geom_point() +
  scale_color_manual(values = color_values) +
  scale_alpha_manual(values = alpha_values) +
  scale_size_manual(values = size_values) +
  scale_shape_manual(values = shape_values) +
  labs(x = paste0("Axis 1 (",round(pcoa_rf$values$Relative_eig[1]*100, 2), " %)"),
       y = paste0("Axis 2 (",round(pcoa_rf$values$Relative_eig[2]*100, 2), " %)")) +
  mytheme

ggsave("forests/chaillou/chaillou-pcoa-rf.png", width = 7.5, height = 5, dpi = "retina")
