library(curatedMetagenomicData)
library(correlationtree)
library(StructFDR)
library(tidyverse)
library(yatah)


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

clades_to_keep <-
  df_abund %>%
  select(-taxonomy) %>%
  gather(key = sample, value = Count, -clade) %>%
  group_by(clade) %>%
  summarise(P = mean(Count > 0), N = sum(Count)) %>%
  arrange(P) %>%
  filter(P > 0.05) %>%
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
  taxtree(lineage_length = mean_lineage_length(tree_cor))

# Random trees 
set.seed(42)
tree_rand_cor <- shuffle_tiplabels(tree_cor)
tree_rand_tax <- shuffle_tiplabels(tree_tax) 


#### Tree FDR ####

## Functions

perm.func <- function (X, Y, ...) {
  return(list(X = X, Y = sample(Y)))
}

test.func.kw <- function (X, Y, O) {
  obj <- apply(X, 1, function(x) {
    kruskal.test(x/O ~ Y)$p.value
  })
  return(list(p.value = obj, e.sign = NULL))
}

my_TreeFDR <- partial(TreeFDR2, B = 200, q.cutoff = 0.5, eff.sign = FALSE,
                      test.func = test.func.kw, perm.func = perm.func)

# set.seed(42)

## Formating

X <-
  df_abund %>%
  select(-taxonomy) %>%
  column_to_rownames("clade") %>%
  as.matrix()

Y <- as_factor(df_sample$study_condition)
O <- df_sample$count

## TreeFDR
fdrobj_cor <- my_TreeFDR(X = X, Y = Y, O = O, tree = tree_cor)
fdrobj_tax <- my_TreeFDR(X = X, Y = Y, O = O, tree = tree_tax)
fdrobj_rand_cor <- my_TreeFDR(X = X, Y = Y, O = O, tree = tree_rand_cor)
fdrobj_rand_tax <- my_TreeFDR(X = X, Y = Y, O = O, tree = tree_rand_tax)

## Agregation

list_fdrobj <- list(fdrobj_cor, fdrobj_tax, fdrobj_rand_cor, fdrobj_rand_tax)

df_pvalues <-
  tibble(genus = names(fdrobj_cor$p.unadj),
         p_raw = fdrobj_cor$p.unadj,
         p_bh = p.adjust(p_raw, method = "BH"),
         p_by = p.adjust(p_raw, method = "BY"),
         p_cor = fdrobj_cor$p.adj,
         p_tax = fdrobj_tax$p.adj,
         p_rand_cor = fdrobj_rand_cor$p.adj,
         p_rand_tax = fdrobj_rand_tax$p.adj) %>%
  gather(-genus, key = "method", value = "pvalue") %>%
  mutate(method = str_remove(method, "p_"),
         method = as_factor(method))

#### Results ####

# Matrix of parameters
matrix(c(map_dbl(list_fdrobj, "k"), map_dbl(list_fdrobj, "rho")),
       ncol = length(list_fdrobj), byrow = FALSE,
       dimnames = list(c("k", "rho"), c("cor", "tax", "rand_cor", "rand_tax")))

# Number of detected genera by method
df_pvalues %>%
  filter(pvalue < 0.05) %>%
  dplyr::count(method)


#### Plot ###

l_pvalues <-
  df_pvalues %>%
  filter(method %in% c("bh", "cor","tax", "rand_cor", "rand_tax")) %>%
  group_by(method) %>%
  group_split() %>%
  map(pull, pvalue) %>%
  set_names(c("bh", "cor", "tax", "rand_cor", "rand_tax"))

color_values <- c("Correlation" = "#C77CFF", "Taxonomy" = "#F8766D",
                  "Random correlation" = "#7CAE00", "Random Taxonomy" = "#FFA500",
                  "BH" = "#4169E1")

linetype_values <- c("Correlation" = "solid", "Taxonomy" = "dotdash",
                     "Random correlation" = "dashed", "Random Taxonomy" = "twodash",
                     "BH" = "dotted")

df_roc <-
  tibble(threshold = seq(0, 0.16, by = 10^-4)) %>%
  mutate(bh = map_dbl(threshold, ~ sum(l_pvalues$bh < .)),
         cor = map_dbl(threshold, ~ sum(l_pvalues$cor < .)),
         tax = map_dbl(threshold, ~ sum(l_pvalues$tax < .)),
         rand_cor = map_dbl(threshold, ~ sum(l_pvalues$rand_cor < .)),
         rand_tax = map_dbl(threshold, ~ sum(l_pvalues$rand_tax < .))) %>%
  gather(key = "method", value = "number", -threshold) %>%
  mutate(method = factor(method, levels = c("bh", "cor", "tax", "rand_cor", "rand_tax"),
                         labels = c("BH", "Correlation", "Taxonomy",
                                    "Random correlation", "Random Taxonomy")),
         method = fct_reorder2(method, threshold, number))


df_arrow <- tribble(~x1, ~y1, ~x2, ~y2,
                    0.063, 13, 0.0508, 13.8,
                    0.037, 17, 0.0492, 16.2)

ggplot(df_roc) +
  aes(x = threshold, y = number) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  geom_line(aes(color = method, linetype = method), size = 1.5, key_glyph = "timeseries") +
  geom_curve(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.1, "inch")),
             size = 0.4, curvature = -0.2, color = "grey20") +
  geom_text(x = 0.032, y = 17, label = 16, size = 7.5) +
  geom_text(x = 0.068, y = 13, label = 14, size = 7.5) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  labs(x = "Threshold", y = "Number of detected genera",
       color = "Method", linetype = "Method") +
  coord_cartesian(clip = 'off') +
  theme_minimal() +
  theme(legend.position = c(0.745, 0.3),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 19))

ggsave("real_datasets/zeller_genus/zeller_genus-number_detected.png", width = 7.5, height = 5, dpi = "retina")

