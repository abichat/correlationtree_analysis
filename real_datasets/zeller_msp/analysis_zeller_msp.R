library(correlationtree)
library(tidyverse)
library(ggstance)
library(cowplot)
library(ggtree)

#### Data ####

## Sample 

df_sample <- read_csv("real_datasets/zeller_msp/condition.csv", col_types = "cc")

## Abundance 

df_abund <- read_tsv("real_datasets/zeller_msp/zeller_all.MSP.sum.txt", skip = 3, 
                     col_types = cols(.default = col_double(), 
                                      MSP = col_character()))
## Filtering 

prevalence_min <- 0.05

df_abund %>%
  gather(key = sample, value = count, -MSP) %>%
  group_by(MSP) %>%
  summarise(P = mean(count > 0), N = sum(count)) %>%
  arrange(P) %>%
  filter(P <= prevalence_min) %>%
  pull(MSP) %>% 
  unique() # All msp are above the minimum prevalence

## Correlation tree

set.seed(42)
tree_cor <- correlation_tree(df_abund, method = "spearman")
tree_rand_cor <- shuffle_tiplabels(tree_cor)


#### Tree FDR ####

## Functions

perm.func <- function (X, Y, ...) {
  return(list(X = X, Y = sample(Y)))
}

test.func.kw <- function (X, Y) {
  obj <- apply(X, 1, function(x) {
    kruskal.test(x ~ Y)$p.value
  })
  return(list(p.value = obj, e.sign = NULL))
}

my_TreeFDR <- partial(TreeFDR2, B = 200, q.cutoff = 0.5, eff.sign = FALSE,
                      test.func = test.func.kw, perm.func = perm.func)

## Formating 

X <-
  df_abund %>%
  column_to_rownames("MSP") %>%
  as.matrix()

X <- X[, df_sample$sample]

Y <- as_factor(df_sample$study_condition)

## Tree FDR

set.seed(42)

fdrobj_cor <- my_TreeFDR(X = X, Y = Y, tree = tree_cor)
fdrobj_rand_cor <- my_TreeFDR(X = X, Y = Y, tree = tree_rand_cor)

## Agregation

list_fdrobj <- list(fdrobj_cor, fdrobj_rand_cor)

df_pvalues <-
  tibble(genus = names(fdrobj_cor$p.unadj),
         p_raw = fdrobj_cor$p.unadj,
         p_bh = p.adjust(p_raw, method = "BH"),
         p_by = p.adjust(p_raw, method = "BY"),
         p_cor = fdrobj_cor$p.adj,
         p_rand_cor = fdrobj_rand_cor$p.adj) %>%
  gather(-genus, key = "method", value = "pvalue") %>%
  mutate(method = str_remove(method, "p_"),
         method = as_factor(method))

#### Results ####

# Matrix of parameters

matrix(c(map_dbl(list_fdrobj, "k"), map_dbl(list_fdrobj, "rho")),
       ncol = length(list_fdrobj), byrow = TRUE,
       dimnames = list(c("k", "rho"), c("cor", "rand_cor")))

# Number of detected genera by method
df_pvalues %>% 
  filter(pvalue < 0.05) %>% 
  count(method)

#### Plots ####

## Number of detected MSP

color_values <- c("Correlation" = "#C77CFF", "Taxonomy" = "#F8766D",
                  "Random correlation" = "#7CAE00", "Random Taxonomy" = "#FFA500",
                  "BH" = "#4169E1")

linetype_values <- c("Correlation" = "solid", "Taxonomy" = "dotdash",
                     "Random correlation" = "dashed", "Random Taxonomy" = "twodash",
                     "BH" = "dotted")

l_pvalues <-
  df_pvalues %>%
  group_by(method) %>%
  group_split() %>%
  map(pull, pvalue) %>%
  set_names(c("raw", "bh", "by", "cor", "rand_cor"))

df_roc <-
  tibble(threshold = seq(0, 0.15, by = 10^-4)) %>%
  mutate(BH = map_dbl(threshold, ~ sum(l_pvalues$bh < .)),
         Correlation = map_dbl(threshold, ~ sum(l_pvalues$cor < .)),
         `Random correlation` = map_dbl(threshold, ~ sum(l_pvalues$rand_cor < .))) %>%
  gather(key = "method", value = "number", -threshold) %>%
  mutate(method = as_factor(method))

df_arrow <- tribble(~x1, ~y1, ~x2, ~y2,
                    0.058, 68, 0.0508, 76.8,
                    0.063, 84, 0.0508, 84.8,
                    0.037, 91, 0.0492, 90.2)

ggplot(df_roc) +
  aes(x = threshold, y = number) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  geom_line(aes(color = method, linetype = method), size = 1.5, key_glyph = "timeseries") +
  geom_curve(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.1, "inch")),
             size = 0.4, curvature = -0.2, color = "grey20") +
  geom_text(x = 0.063, y = 68, label = 77, size = 7.5) +
  geom_text(x = 0.068, y = 84, label = 85, size = 7.5) +
  geom_text(x = 0.032, y = 91, label = 90, size = 7.5) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  labs(x = "Threshold", y = "Number of detected MSPs",
       color = "Method", linetype = "Method") +
  coord_cartesian(clip = 'off') +
  theme_minimal() +
  theme(legend.position = c(0.745, 0.3),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 19))

ggsave("real_datasets/zeller_msp/zeller_msp-number_detected.png", 
       width = 7.5, height = 5, dpi = "retina")

## Difference of p-values

tbl_pvalue <-
  df_pvalues %>% 
  spread(key = method, value = pvalue) %>% 
  mutate(Detected = case_when(cor <= 0.05 & bh <= 0.05 ~ "Correlation and BH",
                              cor <= 0.05              ~ "Correlation",
                              bh <= 0.05               ~ "BH",
                              TRUE                     ~ "None"))

color_detected <- c("Correlation and BH" = "forestgreen", "Correlation" = "blue", "None" = alpha("grey80", 0.4))

p_facet <-
  tree_cor %>% 
  ggtree(color = "grey30") %>% 
  facet_plot(panel = "Difference between correlation tree and BH corrected p-values", data = tbl_pvalue, geom = geom_barh, 
             mapping = aes(x = cor-bh, color = Detected, alpha = Detected), 
             stat = "identity", show.legend = FALSE) %<+%
  tbl_pvalue +
  geom_tippoint(aes(subset = Detected != "None", color = Detected), size = 3) +
  scale_alpha_manual(values = c("None" = 0.2, "Correlation and BH" = 1, "Correlation" = 1)) +
  scale_color_manual(values = color_detected, name = "Detected by") +
  theme(legend.position = "bottom", text = element_text(size = 22))

p_facet$data$.panel <- factor(p_facet$data$.panel, 
                              levels = levels(p_facet$data$.panel),
                              labels = c("Correlation tree", "Difference between correlation tree and BH corrected p-values"))

p_facet_full <-
  p_facet + 
  geom_segment(data = data.frame(.panel = "Difference between correlation tree and BH corrected p-values", 
                                 xintercept = 0, ymin = 0, ymax = 879),
               aes(x = xintercept, xend = xintercept, y = ymin, yend = ymax), color = "grey50")

p_scatter_all <- 
  tbl_pvalue %>% 
  ggplot() + 
  aes(x = bh, y = cor, color = Detected) + 
  geom_rect(xmin = 0, ymin = 0, xmax = 0.06, ymax = 0.06, alpha = 0.5, 
            fill = "tan", color = NA) +
  geom_abline() + 
  geom_vline(xintercept = 0.05) + 
  geom_hline(yintercept = 0.05) +
  geom_point(size = 4) +  
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_color_manual(values = color_detected) +
  labs(x = "BH corrected p-value", y = "Correlation tree corrected p-value") +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 19))

p_scatter_zoom <- 
  tbl_pvalue %>% 
  filter(Detected != "None") %>% 
  ggplot() + 
  aes(x = bh, y = cor, color = Detected) + 
  geom_abline() + 
  geom_vline(xintercept = 0.05) + 
  geom_hline(yintercept = 0.05) +
  geom_point() + 
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 0.056)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 0.056)) +
  scale_color_manual(values = color_detected) +
  labs(x = "BH corrected p-value", y = "Correlation tree corrected p-value") +
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 19))

plot_grid(p_facet_full, 
          plot_grid(p_scatter_all, p_scatter_zoom, ncol = 2), ncol = 1)

ggsave("real_datasets/zeller_msp/zeller_msp-pvalues.png", width = 15, height = 15, dpi = "retina")
