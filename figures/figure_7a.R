library(tidyverse)

source("figures/theme.R")

list_fdrobj <- readRDS("real_datasets/zeller_genus/zeller_genus-fdrobj.rds")

l_pvalues <-
  tibble(genus = names(list_fdrobj$cor$p.unadj),
         p_raw = list_fdrobj$cor$p.unadj,
         p_bh = p.adjust(p_raw, method = "BH"),
         p_by = p.adjust(p_raw, method = "BY"),
         p_cor = list_fdrobj$cor$p.adj,
         p_tax = list_fdrobj$tax$p.adj,
         p_rand_cor = list_fdrobj$rand_cor$p.adj,
         p_rand_tax = list_fdrobj$rand_tax$p.adj) %>%
  gather(-genus, key = "method", value = "pvalue") %>%
  mutate(method = str_remove(method, "p_"),
         method = as_factor(method)) %>%
  filter(method %in% c("bh", "cor","tax", "rand_cor", "rand_tax")) %>%
  group_by(method) %>%
  group_split() %>%
  map(pull, pvalue) %>%
  set_names(c("bh", "cor", "tax", "rand_cor", "rand_tax"))

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
                                    "Random Correlation", "Random Taxonomy")),
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
  geom_text(x = 0.068, y = 13, label = 14, size = 7.5) +
  geom_text(x = 0.032, y = 17, label = 16, size = 7.5) +
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

ggsave("figures/figure_7a.png", width = 7.5, height = 5, dpi = "retina")
