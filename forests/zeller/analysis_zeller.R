library(curatedMetagenomicData)
library(correlationtree)
library(tidyverse)
library(broom)
library(furrr)
library(yatah)
library(ape)

# Parallelize the code
options("future.fork.enable" = TRUE)
plan(multicore)


#### Data ####

exprSet <-
  "ZellerG_2014.metaphlan_bugs_list.stool" %>%
  curatedMetagenomicData(dryrun = FALSE, counts = TRUE) %>%
  mergeData()


## Samples 

df_sample <-
  exprSet %>%
  pData() %>%
  as_tibble(rownames = "sample") %>%
  select(sample, study_condition)

## Abundance & taxonomy 

df_abund <-
  exprSet %>%
  exprs() %>%
  as_tibble(rownames = "taxonomy") %>%
  filter(is_rank(taxonomy, "genus"),
         is_clade(taxonomy, "Bacteria", "kingdom")) %>%
  mutate(clade = last_clade(taxonomy)) %>%
  select(clade, taxonomy, everything())

total_count <-
  df_abund %>%
  select(-clade, -taxonomy) %>%
  colSums()

df_sample$count <- total_count

prevalence_min <- 0.05

clades_to_keep <-
  df_abund %>%
  select(-taxonomy) %>%
  gather(key = sample, value = Count, -clade) %>%
  group_by(clade) %>%
  summarise(P = mean(Count > 0), N = sum(Count)) %>%
  arrange(P) %>%
  filter(P > prevalence_min) %>%
  pull(clade)

df_abund <- filter(df_abund, clade %in% clades_to_keep)

## Trees

# Correlation tree
tree_cor <-
  df_abund %>%
  select(-taxonomy) %>%
  correlation_tree(method = "spearman")

# Taxonomy
tree_tax <-
  df_abund %>%
  pull(taxonomy) %>%
  taxtable() %>%
  taxtree(lineage_length = mean_lineage_length(tree_cor)) %>% 
  multi2di()

tree_tax$node.label <- NULL


#### Forest ####

set.seed(42)

N_boot <- 250 # To ensure about 100 bootstrapped trees
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

tree_types <- 
  factor(c("Correlation", "Taxonomy", 
           rep("Bootstrap", N_boot), 
           rep("Random Correlation", N_rand), 
           rep("Random Taxonomy", N_rand)),
         levels = c("Correlation", "Bootstrap", 
                    "Random Correlation", "Random Taxonomy", "Taxonomy"))


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


color_values <- c("Correlation" = "#C77CFF", "Taxonomy" = "#F8766D", 
                  "Bootstrap" = "#00BFC4", "Random Correlation" = "#7CAE00",
                  "Random Taxonomy" = "#FFA500")

size_values <- c("Correlation" = 4, "Taxonomy" = 4, "Bootstrap" = 1, 
                 "Random Correlation" = 1, "Random Taxonomy" = 1)

alpha_values <- c("Correlation" = 0.8, "Taxonomy" = 0.8, "Bootstrap" = .4, 
                  "Random Correlation" = .4, "Random Taxonomy" = .4)

shape_values <- c("Correlation" = 17, "Taxonomy" = 16, "Bootstrap" = 2, 
                  "Random Correlation" = 6, "Random Taxonomy" = 1)


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

ggsave("forests/zeller/zeller-bhv-boxplot.png", width = 7.5, height = 5, dpi = "retina")

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

ggsave("forests/zeller/zeller-rf-boxplot.png", width = 7.5, height = 5, dpi = "retina")


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

ggsave("forests/zeller/zeller-bhv-pcoa.png", width = 7.5, height = 5, dpi = "retina")

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

ggsave("forests/zeller/zeller-rf-pcoa.png", width = 7.5, height = 5, dpi = "retina")
