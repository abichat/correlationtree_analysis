library(correlationtree)
library(tidyverse)

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

beepr::beep(5)

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

#### Plot ####

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

# df_roc_005 <- 
#   df_roc %>% 
#   filter(threshold == 0.05) %>% 
#   mutate(position = number + c(-5, 5))
# 
# ggplot(df_roc) +
#   aes(x = threshold, y = number) +
#   geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
#   geom_line(aes(color = method, linetype = method), size = 1.5, key_glyph = "timeseries") +
#   geom_label(data = df_roc_005, aes(color = method, label = number, y = position), 
#              size = 7.5, show.legend = FALSE) +
#   scale_color_manual(values = color_values) +
#   scale_linetype_manual(values = linetype_values) +
#   labs(x = "Threshold", y = "Number of detected MSPs",
#        color = "Method", linetype = "Method") +
#   coord_cartesian(clip = 'off') +
#   theme_minimal() +
#   theme(legend.position = c(0.745, 0.3),
#         axis.title = element_text(size = 22),
#         axis.text = element_text(size = 15),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 19))

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
