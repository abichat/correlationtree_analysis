library(correlationtree)
library(tidyverse)
library(structSSI)
library(igraph)
library(ape)


#### Data ####

data("chlamydiae")

abundances <- as.data.frame(chlamydiae@otu_table)
OTU <- rownames(abundances)
environments <- sample_data(chlamydiae)$SampleType

tree_phy <- chlamydiae@phy_tree
tree_cor <- correlation_tree(abundances, col = 0, method = "spearman", remove = FALSE)


#### Hierarchical FDR ####

alpha <- 0.1

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
head(EL_cor)

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

## Results

# A posteriori FDR control for phylogenetic correction
1.44 * EstimatedHFDRControl(hpval_phy)$tip
# A posteriori FDR control for correlation correction
1.44 * EstimatedHFDRControl(hpval_cor)$tip

# Number of detected species by phylogenetic correction
length(detected_phy)
# Number of detected species by correlation correction
length(detected_cor)
# Number of detected species by BH with similar FDR
tbl_pvalues %>% filter(bh < 0.324) %>% pull(OTU) %>% length()


#### Plot ####

library(cowplot)
library(ggtree)
library(scales)

## Trees 

p_phy_phy <-
  tree_phy %>%
  ggtree(branch.length="none", color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(phy)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = phy < 0.1), vjust = -0.3, hjust = 1, 
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1) 

p_phy_cor <-
  tree_phy %>% 
  ggtree(branch.length="none", color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(cor)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = cor < 0.1), vjust = -0.3,
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1) +
  scale_x_reverse()

p_cor_phy <-
  tree_cor %>% 
  ggtree(color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(phy)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = phy < 0.1), vjust = -0.3, hjust = 1,
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1)

p_cor_cor <-
  tree_cor %>% 
  ggtree(color = "grey30") %<+%
  tbl_pvalues +
  geom_tippoint(aes(size = -log10(cor)), color = "grey60", alpha = 0.5) +
  geom_tiplab(aes(color = cor < 0.1), vjust = -0.3,
              size = 2.5, fontface = "bold") +
  scale_color_viridis_d(direction = -1) +
  scale_x_reverse()


p_trees <-
  plot_grid(p_phy_phy, p_phy_cor, NULL, NULL, p_cor_phy, p_cor_cor, 
            ncol = 2, rel_heights = c(1, 0.1, 1),
            labels = c("A", "B", "", "", "C", "D"), 
            label_x = c(0, 0.9, 0, 0, 0, 0.9), label_y = c(1, 1, 0, 0, 0.12, 0.12))

## Boxplots 

otu_cor <-
  tbl_pvalues %>% 
  filter(cor < 0.1, phy > 0.1 | is.na(phy)) %>% 
  pull(OTU)

df_abund_cor <-
  abundances %>% 
  as.data.frame() %>% 
  as_tibble(rownames = "OTU") %>% 
  filter(OTU %in% otu_cor) %>% 
  gather(key = "X.SampleID", value = "Abundance", -OTU) %>%
  left_join(tibble("X.SampleID" = colnames(abundances),  SampleType = environments), 
            by = "X.SampleID") %>% 
  mutate(SampleType = fct_recode(SampleType, FE = "Feces", FW = "Freshwater", 
                                 CK ="Freshwater (creek)", MO = "Mock", 
                                 OC = "Ocean", SE = "Sediment (estuary)", 
                                 SK = "Skin", SO = "Soil", TO = "Tongue"),
         SampleType = fct_relevel(SampleType, "SO", "SE", "OC", "CK", "FW",
                                  "SK", "TO", "FE", "MO"),
         OTU = as_factor(OTU)) 

mytheme <-
  theme_minimal() +
  theme(strip.text = element_text(size = 10),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) 

p_abund1 <- 
  df_abund_cor %>% 
  filter(OTU == 547579) %>% 
  ggplot() +
  aes(x = SampleType, y = Abundance) +
  geom_boxplot() +
  facet_wrap(~ OTU, scales = "free_y") +
  scale_y_sqrt(breaks = pretty_breaks(4)) +
  labs(x = NULL) +
  mytheme

p_abund2 <- 
  df_abund_cor %>% 
  filter(OTU != 547579) %>% 
  ggplot() +
  aes(x = SampleType, y = Abundance) +
  geom_boxplot() +
  facet_wrap(~ OTU, scales = "free_y") +
  scale_y_sqrt(breaks = pretty_breaks(4)) +
  labs(x = NULL, y = NULL) +
  mytheme

## Full plot 

p_trees <-
  plot_grid(p_phy_phy, p_phy_cor, NULL, NULL, p_cor_phy, p_cor_cor, 
            ncol = 2, rel_heights = c(1, 0.07, 1),
            labels = c("A", "B", "", "", "C", "D"), label_size = 14,
            label_x = c(0, 0.9, 0, 0, 0, 0.9), label_y = c(1, 1, 0, 0, 0.12, 0.12))

ggdraw(p_trees) +
  draw_plot(p_abund1, x = 0.01, y = 0.32, width = 0.3, height = 0.27) +
  draw_plot(p_abund2, x = 0.69, y = 0.32, width = 0.3, height = 0.27) +
  draw_plot_label(label = c("E", "F"), size = 14,
                  x = c(0.3, 0.66), y = c(0.55, 0.55)) +
  draw_plot_label("*", x = 0.59, y = c(0.442, 0.676), colour = "red")

ggsave("real_datasets/chlamydiae/chlamydiae-evidences_on_trees.png", width = 15, height = 12, units = "cm")
