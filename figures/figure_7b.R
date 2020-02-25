library(tidyverse)

source("figures/theme.R")

list_fdrobj <- readRDS("real_datasets/zeller_msp/zeller_msp-fdrobj.rds")

l_pvalues <-
  tibble(genus = names(list_fdrobj$cor$p.unadj),
         p_raw = list_fdrobj$cor$p.unadj,
         p_bh = p.adjust(p_raw, method = "BH"),
         p_by = p.adjust(p_raw, method = "BY"),
         p_cor = list_fdrobj$cor$p.adj,
         p_rand_cor = list_fdrobj$rand_cor$p.adj) %>%
  gather(-genus, key = "method", value = "pvalue") %>%
  mutate(method = str_remove(method, "p_"),
         method = as_factor(method)) %>% 
  group_by(method) %>%
  group_split() %>%
  map(pull, pvalue) %>%
  set_names(c("raw", "bh", "by", "cor", "rand_cor"))

df_roc <-
  tibble(threshold = seq(0, 0.15, by = 10^-4)) %>%
  mutate(BH = map_dbl(threshold, ~ sum(l_pvalues$bh < .)),
         Correlation = map_dbl(threshold, ~ sum(l_pvalues$cor < .)),
         `Random Correlation` = map_dbl(threshold, ~ sum(l_pvalues$rand_cor < .))) %>%
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

ggsave("figures/figure_7b.png", width = 7.5, height = 5, dpi = "retina")
