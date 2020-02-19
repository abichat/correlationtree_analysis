library(correlationtree)
library(biomformat)
library(tidyverse)
library(broom)
library(furrr)
library(yatah)
library(ape)

# Parallelize the code
options("future.fork.enable" = TRUE)
plan(multiprocess)


#### Data ####

biom <- read_biom("forests/ravel/ravel.biom")

## Taxonomy

taxtable <-
  biom %>% 
  observation_metadata() %>% 
  as("matrix") %>% 
  as_tibble(rownames = "OTU") %>% 
  filter(taxonomy1 %in% "Bacteria") %>% 
  select(OTU, kingdom = taxonomy1, phylum = taxonomy2, class = taxonomy3, 
         order = taxonomy4, family = taxonomy5, genus = taxonomy6) %>% 
  filter(genus != "NA") %>% 
  unique()

## Abundance and filtering

prevalence_min <- 0.05

df_abund <- 
  biom %>% 
  biom_data() %>% 
  as("matrix") %>% 
  as_tibble(rownames = "OTU")

df_abund <- 
  taxtable %>% 
  left_join(df_abund, by = "OTU") %>%
  select(-c(OTU, kingdom, phylum, class, order, family)) %>% 
  group_by(genus) %>% 
  summarise_all(sum)

df_abund <-
  df_abund %>% 
  gather(key = "sample", value = "abundance", -genus) %>% 
  group_by(genus) %>% 
  mutate(P = mean(abundance > 0)) %>% 
  filter(P > prevalence_min) %>% 
  select(-P) %>% 
  ungroup() %>% 
  spread(key = sample, value = abundance)

taxtable <- filter(taxtable, genus %in% df_abund$genus)

## Trees

tree_cor <- correlation_tree(df_abund, method = "spearman")
mean_lineage <- mean_lineage_length(tree_cor)

tree_tax <- 
  taxtable %>% 
  select(-OTU) %>% 
  unique() %>% 
  taxtree(lineage_length = mean_lineage) %>% 
  multi2di()
tree_tax$node.label <- NULL


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

trees_rand_tax <-
  N_rand %>%
  rerun(shuffle_tiplabels(tree_tax)) %>%
  map(multi2di) %>%
  reduce(c)

## Aggregation

forest <- c(tree_cor, tree_tax, trees_boot, trees_rand_cor, trees_rand_tax)

tree_labels <- 
  factor(c("Correlation", "Taxonomy", 
           rep("Bootstrap", N_boot), 
           rep("Random Correlation", N_rand), 
           rep("Random Taxonomy", N_rand)),
         levels = c("Correlation", "Bootstrap", 
                    "Random Correlation", "Random Taxonomy", "Taxonomy"))

# saveRDS(tree_labels, "forests/ravel/ravel-tree-labels.rds")
# tree_labels <- readRDS("forests/ravel/ravel-tree-labels.rds")


#### Distances and models ####

# Pairwise distances, could take time
dist_bhv <- future_dist_BHV(forest) # Billera-Holmes-Vogtmann
dist_rf <- dist.topo(unroot(forest)) # Robinson-Foulds

# saveRDS(dist_bhv, "forests/ravel/ravel-dist-bhv.rds")
# saveRDS(dist_rf,  "forests/ravel/ravel-dist-rf.rds")
# dist_bhv <- readRDS("forests/ravel/ravel-dist-bhv.rds")
# dist_rf <-  readRDS("forests/ravel/ravel-dist-rf.rds")

# PCoA
pcoa_bhv <- pcoa(dist_bhv) 
pcoa_rf <- pcoa(dist_rf) 

# saveRDS(pcoa_bhv, "forests/ravel/ravel-pcoa-bhv.rds")
# saveRDS(pcoa_rf,  "forests/ravel/ravel-pcoa-rf.rds")
# pcoa_bhv <- readRDS("forests/ravel/ravel-pcoa-bhv.rds")
# pcoa_rf <-  readRDS("forests/ravel/ravel-pcoa-rf.rds")

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

## Theme

source("figures/theme.R")

## Boxplots

dist_tax_bhv <- filter(dist_df_bhv, Type == "Taxonomy")$Distance
dist_tax_rf  <- filter(dist_df_rf, Type == "Taxonomy")$Distance

dist_df_bhv %>% 
  filter(!Type %in% c("Correlation", "Taxonomy")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_tax_bhv, color = color_values["Taxonomy"]) +
  geom_text(x = 1, y = dist_tax_bhv, vjust = 1.2, size = 4,
            color = color_values["Taxonomy"], label = "Taxonomy") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme

ggsave("forests/ravel/ravel-boxplot-bhv.png", width = 7.5, height = 5, dpi = "retina")

dist_df_rf %>% 
  filter(!Type %in% c("Correlation", "Taxonomy")) %>% 
  ggplot() +
  aes(x = Type, y = Distance) +
  geom_violin(aes(fill = Type), alpha = 0.5) +
  geom_boxplot(aes(color = Type), alpha = 0, notch = TRUE) +
  geom_hline(yintercept = dist_tax_rf, color = color_values["Taxonomy"]) +
  geom_text(x = 1, y = dist_tax_rf, vjust = 1.2, size = 4,
            color = color_values["Taxonomy"], label = "Taxonomy") +
  geom_hline(yintercept = 0, color = color_values["Correlation"]) +
  geom_hline(yintercept = 0, alpha = 0) +
  scale_fill_manual(values = color_values) +
  scale_color_manual(values = color_values) +
  labs(x = NULL, y = "Distance to correlation tree") +
  mytheme

ggsave("forests/ravel/ravel-boxplot-rf.png", width = 7.5, height = 5, dpi = "retina")


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

ggsave("forests/ravel/ravel-pcoa-bhv.png", width = 7.5, height = 5, dpi = "retina")

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

ggsave("forests/ravel/ravel-pcoa-rf.png", width = 7.5, height = 5, dpi = "retina")
