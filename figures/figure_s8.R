library(tidyverse)
library(ggstance)
library(cowplot)
library(ggtree)

tree_cor <- readRDS("real_datasets/zeller_msp/zeller_msp-tree_cor.rds")
list_fdrobj <- readRDS("real_datasets/zeller_msp/zeller_msp-fdrobj.rds")

df_pvalues <-
  tibble(genus = names(list_fdrobj$cor$p.unadj),
         p_raw = list_fdrobj$cor$p.unadj,
         p_bh = p.adjust(p_raw, method = "BH"),
         p_by = p.adjust(p_raw, method = "BY"),
         p_cor = list_fdrobj$cor$p.adj,
         p_rand_cor = list_fdrobj$rand_cor$p.adj) %>%
  gather(-genus, key = "method", value = "pvalue") %>%
  mutate(method = str_remove(method, "p_"),
         method = as_factor(method))

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

ggsave("figures/figure_s8.png", width = 15, height = 15, dpi = "retina")
