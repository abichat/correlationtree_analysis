library(correlationtree)
library(biomformat)
library(structSSI)
library(tidyverse)
library(phyloseq)
library(cowplot)
library(ggtree)
library(igraph)
library(broom)
library(ape)


#### Data ####

biom <- read_biom("real_datasets/chaillou/chaillou.biom")

## Taxonomy 

genus_peptoniphilus <- c("Anaerococcus", "Finegoldia", "Gallicola", "Peptoniphilus")
genus_tissierellaceae <- "Tissierella"

taxtable <-
  biom %>% 
  observation_metadata() %>% 
  as("matrix") %>% 
  as_tibble(rownames = "OTU") %>% 
  mutate_at(vars(matches("taxonomy*")), ~ str_sub(., start = 4)) %>% 
  mutate_at(vars(matches("taxonomy*")), ~ str_replace_all(., " ", "_")) %>% 
  mutate_at(vars(matches("taxonomy*")), ~ str_replace_all(., "\\.", "_")) %>% 
  select(OTU, kingdom = taxonomy1, phylum = taxonomy2, class = taxonomy3, 
         order = taxonomy4, family = taxonomy5, genus = taxonomy6) %>% 
  # Fix taxonomy errors
  mutate(family = case_when(genus %in% genus_peptoniphilus   ~ "Peptoniphilaceae",
                            genus %in% genus_tissierellaceae ~ "Tissierellaceae",
                            TRUE                             ~ family),
         order = case_when(genus %in% genus_peptoniphilus    ~ "Tissierellales",
                           genus %in% genus_tissierellaceae  ~ "Tissierellales",
                           TRUE                             ~ order),
         class = case_when(genus %in% genus_peptoniphilus    ~ "Tissierellia",
                           genus %in% genus_tissierellaceae  ~ "Tissierellia",
                           TRUE                             ~ class)) %>% 
  filter(kingdom %in% "Bacteria") %>% 
  select(-genus) %>% 
  filter(family != "NA") %>% 
  unique()

## Filter abundances

prevalence_min <- 0.05

abundances <- 
  biom %>% 
  biom_data() %>% 
  as("matrix") %>% 
  as_tibble(rownames = "OTU")

abundances <-
  abundances %>%
  left_join(select(taxtable, OTU, phylum), by = "OTU") %>% 
  # Keeping only Bacteroidetes
  filter(phylum == "Bacteroidetes") %>% 
  select(-phylum) %>% 
  gather(key = "sample", value = "abundance", -OTU) %>% 
  group_by(OTU) %>% 
  mutate(P = mean(abundance > 0)) %>% 
  filter(P > prevalence_min) %>% 
  select(-P) %>% 
  ungroup() %>% 
  spread(key = sample, value = abundance) %>%
  as.data.frame() %>% 
  column_to_rownames("OTU") 

OTU <- rownames(abundances)


## Environment

environments <- 
  abundances %>% 
  colnames() %>% 
  sort() %>% 
  str_remove_all("\\..*$") %>% 
  as_factor()


## Trees 

tree_cor <- correlation_tree(abundances, matrix = TRUE, method = "spearman")
mean_lineage <- mean_lineage_length(tree_cor)

# write.tree(tree_cor, file = "real_datasets/chaillou/cortree_chaillou.nwk")

tree_phy <- 
  read.tree("real_datasets/chaillou/phytree_chaillou.nwk") %>% 
  prune_taxa(OTU, .)

tree_phy$edge.length <- tree_phy$edge.length * mean_lineage / mean_lineage_length(tree_phy)


#### Hierarchical FDR ####

alpha <- 0.01

## Phylogeny correction

EL_phy <- 
  tree_phy %>% 
  as.igraph() %>% 
  get.edgelist()

pval_phy <- treePValues(EL_phy, abundances, environments)
hpval_phy <- hFDR.adjust(pval_phy, EL_phy, alpha = alpha)
tbl_phy <- tibble(OTU = OTU, phy = hpval_phy@p.vals[OTU, ]$adjp)

## Correlation correction

EL_cor <- 
  tree_cor %>% 
  as.igraph() %>% 
  get.edgelist()

pval_cor <- treePValues(EL_cor, abundances, environments)
hpval_cor <- hFDR.adjust(pval_cor, EL_cor, alpha = alpha)
tbl_cor <- tibble(OTU = OTU, cor = hpval_cor@p.vals[OTU, ]$adjp)

## Benjamini-Hochberg correction

tbl_notree <- 
  pval_phy[OTU] %>% 
  tibble(OTU = names(.), unadj = .) %>% 
  mutate(bh = p.adjust(unadj, method = "fdr"))

## Aggregation

tbl_pvalues <- reduce(list(tbl_notree, tbl_phy, tbl_cor), left_join, by = "OTU")

detected_phy <-
  tbl_pvalues %>%
  filter(phy < alpha) %>%
  pull(OTU)

detected_cor <-
  tbl_pvalues %>%
  filter(cor < alpha) %>%
  pull(OTU)


#### Results ####

# A posteriori FDR control for phylogenetic correction
1.44 * EstimatedHFDRControl(hpval_phy)$tip
# A posteriori FDR control for correlation correction
1.44 * EstimatedHFDRControl(hpval_cor)$tip

# Number of species detected by phylogenetic correction
length(detected_phy)
# Number of species detected by correlation correction
length(detected_cor)
# Number of species detected by both corrections
length(intersect(detected_cor, detected_phy))
# Number of detected species by BH with similar FDR
tbl_pvalues %>% filter(bh < 0.04) %>% nrow()


#### Plots ####

## Abundances only one

tbl_pvalues_filtered <- 
  tbl_pvalues %>% 
  mutate(Detected = case_when(cor <= 0.01 & phy <= 0.01 ~ "Both",
                              cor <= 0.01               ~ "Correlation",
                              phy <= 0.01               ~ "Phylogeny",
                              TRUE                      ~ "None"),
         Detected = fct_relevel(Detected, "Correlation", "Both", "Phylogeny", "None"),
         OTU_displayed = if_else(Detected != "None", OTU, "", missing = ""))

df_4_boxplots <-
  abundances %>% 
  as_tibble(rownames = "OTU") %>% 
  gather(key = "Sample", value = "Count", -OTU) %>% 
  left_join(select(tbl_pvalues_filtered, OTU, Detected), by = "OTU") %>% 
  arrange(Detected, OTU) %>% 
  mutate(OTU = as_factor(OTU),
         Env = str_sub(Sample, end = 2),
         Env = case_when(
           Env == "BH" ~ "GB",
           Env == "CD" ~ "Sh",
           Env == "DL" ~ "SB",
           Env == "FC" ~ "CF",
           Env == "FS" ~ "SF",
           Env == "MV" ~ "PS",
           Env == "SF" ~ "SS",
           Env == "VH" ~ "GV"),
         Env = fct_relevel(Env, "SB", "PS", "GB", "GV", "Sh", "CF", "SF", "SS"))

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

ggsave("real_datasets/chaillou/chaillou-abund_only_one.png", width = 15, height = 12, units = "cm")

## Evidences on tree

palette_detected <- setNames(c("blue", "purple", "firebrick", "black"), 
                             c("Correlation", "Both", "Phylogeny", "None"))

detected_only_one <- 
  tbl_pvalues_filtered %>% 
  filter(Detected %in% c("Correlation", "Phylogeny")) %>% 
  pull(OTU)

p1 <-
  tree_cor %>%
  ggtree(branch.length="none", color = "grey30") %<+%
  tbl_pvalues_filtered +
  geom_hilight(node = 163, fill = "green", alpha = 0.5) +
  geom_tippoint(aes(size = -log10(cor), color = Detected), alpha = 0.5) +
  geom_tiplab(aes(label = OTU_displayed, color = Detected), vjust = 0, hjust = 1, 
              size = 2.5, fontface = "bold", show.legend = FALSE) +
  guides(size = FALSE) +
  scale_color_manual(values = palette_detected, name = "Detected by")

legend_tree <- get_legend(p1 + theme(legend.position = "bottom"))

p2 <-
  tree_phy %>% 
  ggtree(color = "grey30") %<+%
  tbl_pvalues_filtered +
  geom_tippoint(aes(size = -log10(phy), color = Detected), alpha = 0.5) +
  geom_tiplab(aes(label = OTU_displayed, color = Detected), vjust = 0, hjust = 0,
              size = 2.5, fontface = "bold") +
  scale_color_manual(values = palette_detected) + 
  theme(legend.position = "none") + 
  scale_x_reverse()

plot_grid(
  plot_grid(p1, p2, 
            ncol = 2, 
            labels = c("Correlation tree", "Phylogeny"), 
            label_x = c(0, 0.3), 
            label_y = c(1, 1)),
  legend_tree, 
  ncol = 1, rel_heights = c(2, 0.06))

ggsave("real_datasets/chaillou/chaillou-evidences_on_trees.png", width = 15, height = 22, units = "cm")

## Focus on subtree

otu_subtree <- paste0("otu_0", c("0241", "0516", "0519", "0656", "1495"))
subtree_cor <- keep.tip(tree_cor, otu_subtree)

subtree_cor$node.label <- as.character(6:9)

p_value_lm <- function(otus) {
  glance(lm(colSums(abundances[otu_subtree[otus], ]) ~ environments))$p.value
}

p_value_kw <- function(otus) {
  glance(kruskal.test(colSums(abundances[otu_subtree[otus], ]) ~ environments))$p.value
}

tbl_node_pv <-
  tibble(mrca = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
         merge = list(1, 2, 3, 4, 5, 1:5, 1:2, 3:5, 3:4),
         p_lm = map_dbl(merge, p_value_lm),
         p_kw = map_dbl(merge, p_value_kw)) %>% 
  mutate(mrca = if_else(mrca < 6, otu_subtree[mrca], as.character(mrca)),
         detected = if_else(mrca %in% detected_phy, "Phylogeny", "None")) %>% 
  mutate_at(vars(starts_with("p_")), ~ prettyNum(., format = "e", digits = 2))

tbl_node_pv$p_lm[6] <- paste(tbl_node_pv$p_lm[6], "(F-test p-value)")
tbl_node_pv$p_kw[6] <- paste(tbl_node_pv$p_kw[6], "(KW-test p-value)")

n_y <- 0.08
n_x <- 0.005
hj <- 0
col_lm <- "black"
col_kw <- "grey60"

p_subtree <-
  ggtree(subtree_cor) %<+% 
  tbl_node_pv +
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
  filter(OTU %in% otu_subtree) %>% 
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

ggsave("real_datasets/chaillou/chaillou-focus_subtree.png", width = 15, height = 15, units = "cm")
