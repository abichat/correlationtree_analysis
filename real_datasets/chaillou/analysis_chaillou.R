library(correlationtree)
library(biomformat)
library(structSSI)
library(tidyverse)
library(phyloseq)
library(igraph)
library(ape)
library(ggbeeswarm)
library(cowplot)
library(ggtree)
library(scales)
# setwd(here::here())
source("R_scripts/helpers.R")

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

write.tree(tree_cor, file = "real_datasets/chaillou/cortree_chaillou.nwk")

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
tbl_pvalues %>% filter(bh < 0.04) %>% pull(OTU) %>% length()

## Highlight the 6 taxa that are specific to phy 
otu_phy <- filter(tbl_pvalues, phy <= alpha, cor >= alpha | is.na(cor)) %>% pull(OTU)
otu_cor <- filter(tbl_pvalues, cor <= alpha, phy >= alpha | is.na(phy)) %>% pull(OTU)
otu_test <- paste0("otu_0", c("0241", "0516", "1495", "0519", "0656"))


my_plot_chaillou <- partial(my_plot, data = abundances, env = food_type)

my_plot_chaillou(otu_phy) + ggtitle("OTUs only detected using the phylogeny")
my_plot_chaillou(otu_cor) + ggtitle("OTUs only detected using the correlation")
my_plot_chaillou(otu_test) + ggtitle("OTUs in clade of 01495")  

## Show OTUs on taxonomy
p <- facing_trees(tree_cor, tree_phy, tbl_pvalues, cor, phy, p_thresh = 0.01)

plot(p)
