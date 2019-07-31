library(correlationtree)
library(tidyverse)
library(structSSI)
library(broom)
library(furrr)
library(ape)

# Parallelize the code
options("future.fork.enable" = TRUE)
plan(multiprocess)


#### Data ####

data("chlamydiae")

## Abundance

df_abund <- 
  chlamydiae@otu_table %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "OTU")

## Trees

# Correlation tree
tree_cor <- correlation_tree(df_abund, method = "spearman", remove = FALSE)

# Phylogeny
tree_phy <- chlamydiae@phy_tree
tree_phy$edge.length <- 
  mean_lineage_length(tree_cor) * tree_phy$edge.length / 
  mean_lineage_length(tree_phy)
tree_phy$node.label <- NULL


#### Forest ####

set.seed(42)

N_boot <- 100 
N_rand <- 100

## Bootstraps

correlation_tree <- possibly(correlation_tree, NULL)

trees_boot <- 
  N_boot %>% 
  rerun(sample_boot(df_abund)) %>% 
  future_map(correlation_tree, method = "spearman", fill = TRUE) %>%
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

tree_types <- 
  factor(c("Correlation", "Phylogeny", 
           rep("Bootstrap", N_boot), 
           rep("Random Correlation", N_rand), 
           rep("Random Phylogeny", N_rand)),
         levels = c("Correlation", "Bootstrap", 
                    "Random Correlation", "Random Phylogeny", "Phylogeny"))


#### Distances and models ####

# Pairwise distances, could take time
dist_bhv <- future_dist_BHV(forest) # Billera-Holmes-Vogtmann
dist_rf <- dist.topo(unroot(forest)) # Robinson-Foulds

# PCoA
pcoa_bhv <- pcoa(dist_bhv) 
pcoa_rf <- pcoa(dist_rf) 

# Distances to correlation tree
dist_df_bhv <- tibble(Distance = dist_bhv[1:(length(tree_types) - 1)], Type = tree_types[-1])
dist_df_rf  <- tibble(Distance = dist_rf[1:(length(tree_types) - 1)],  Type = tree_types[-1])

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

mytheme <-
  theme_minimal() +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 12),
        legend.position = "none")


color_values <- c("Correlation" = "#C77CFF", "Phylogeny" = "#F8766D", 
                  "Bootstrap" = "#00BFC4", "Random Correlation" = "#7CAE00",
                  "Random Phylogeny" = "#FFA500")

size_values <- c("Correlation" = 4, "Phylogeny" = 4, "Bootstrap" = 1, 
                 "Random Correlation" = 1, "Random Phylogeny" = 1)

alpha_values <- c("Correlation" = 0.8, "Phylogeny" = 0.8, "Bootstrap" = .4, 
                  "Random Correlation" = .4, "Random Phylogeny" = .4)

shape_values <- c("Correlation" = 17, "Phylogeny" = 16, "Bootstrap" = 2, 
                  "Random Correlation" = 6, "Random Phylogeny" = 1)


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

ggsave("forests/chlamydiae/chlamydiae-bhv-boxplot.png", width = 7.5, height = 5, dpi = "retina")

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

ggsave("forests/chlamydiae/chlamydiae-rf-boxplot.png", width = 7.5, height = 5, dpi = "retina")


## PCoAs

pcoa_bhv$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_types) %>% 
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

ggsave("forests/chlamydiae/chlamydiae-bhv-pcoa.png", width = 7.5, height = 5, dpi = "retina")

pcoa_rf$vectors %>%
  as_tibble() %>%
  mutate(Type = tree_types) %>% 
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

ggsave("forests/chlamydiae/chlamydiae-rf-pcoa.png", width = 7.5, height = 5, dpi = "retina")
