library(tidyverse)
library(cowplot)
library(evabic)

source("figures/theme.R")

Ntaxa <- readRDS("simulations/non_parametric/simus_np-ntaxa.rds")
df_gathered <- readRDS("simulations/non_parametric/simus_np-df_gathered.rds")

df_bh <- 
  df_gathered %>% 
  drop_na(fdr_obj) %>% 
  group_by(time, ID) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  mutate(method = "bh",
         praw = map(fdr_obj, "p.unadj"),
         pbh = map(praw, p.adjust, method = "BH"),
         detected = map(pbh, ~ names(.)[. < 0.05])) %>% 
  mutate(k = NA, rho = NA, smoothing_mean = NA) %>%
  select(fc, nH1, taxa_diffs, method, B, detected, k, rho, smoothing_mean)

df_treefdr <- 
  df_gathered %>% 
  drop_na(fdr_obj) %>% 
  mutate(detected = map(fdr_obj, ~ names(.$p.unadj)[.$p.adj < 0.05])) %>% 
  mutate(k = map_dbl(fdr_obj, "k"), 
         rho = map_dbl(fdr_obj, "rho")) %>% 
  mutate(smoothing_mean = map_dbl(fdr_obj, ~ mean(abs(.$z.adj - .$z.unadj)))) %>% 
  select(fc, nH1, taxa_diffs, method, B, detected, k, rho, smoothing_mean)

df_eval <-
  rbind(df_treefdr, df_bh) %>% 
  mutate(pi0 = 100 * (Ntaxa - nH1) / Ntaxa, 
         tidyebc = map2(detected, taxa_diffs, ebc_tidy, m = Ntaxa,
                        measures = c("BACC", "ACC", "TPR", "FDR", "F1"))) %>% 
  unnest(tidyebc) %>% 
  mutate(BACC = ifelse(is.nan(BACC), 0, BACC)) %>% 
  mutate(FDR = ifelse(is.nan(FDR), 0, FDR)) %>% 
  select(-taxa_diffs, -detected)

df_ebc <- 
  df_eval %>% 
  arrange(nH1) %>% 
  mutate(nH1 = as_factor(nH1)) %>% 
  mutate(method = factor(method, levels = c("bh", "cor", "tax", "randcor", "randtax"), 
                         labels = c("BH", "Correlation", "Taxonomy",
                                    "Random Correlation", "Random Taxonomy")))

df_TPR <-
  df_ebc %>% 
  group_by(pi0, fc, method) %>% 
  summarise(mean = mean(TPR), sd = sd(TPR), count = n()) %>% 
  arrange(desc(method)) %>% 
  mutate(infbound = mean - sd/sqrt(count), 
         supbound = mean + sd/sqrt(count))

labels <- c("5" = "fc = 5", "10" = "fc = 10", "15" = "fc = 15", "20" = "fc = 20")

p_TPR <-
  ggplot(df_TPR) +
  aes(x = pi0, y = mean, color = method) +
  geom_errorbar(aes(ymin = infbound, ymax = supbound, width = 1)) +
  geom_line(aes(linetype = method, group = method)) +
  geom_point(show.legend = FALSE) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  facet_wrap(~ fc, ncol = 4, labeller = labeller(fc = labels)) +
  labs(x = NULL, y = "TPR",
       color = "Method", linetype = "Method") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8), 
        legend.background = element_rect("transparent"))

legend <- get_legend(p_TPR)

p_TPR <- p_TPR + theme(legend.position = "none")

df_FDR <-
  df_ebc %>% 
  group_by(pi0, fc, method) %>% 
  summarise(mean = mean(FDR), sd = sd(FDR), count = n()) %>% 
  arrange(desc(method)) %>% 
  mutate(infbound = mean - sd/sqrt(count), 
         supbound = mean + sd/sqrt(count))

p_FDR <-
  ggplot(df_FDR) +
  aes(x = pi0, y = mean, color = method) +
  geom_hline(yintercept = 0.05, color = "red", alpha = 0.6) +
  geom_errorbar(aes(ymin = infbound, ymax = supbound, width = 1)) +
  geom_line(aes(linetype = method, group = method)) +
  geom_point(show.legend = FALSE) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  facet_wrap(~ fc, ncol = 4, labeller = labeller(fc = labels)) +
  labs(x = "Proportion of null hypothesis", y = "FDR",
       color = "Method", linetype = "Method") +
  theme_bw() +
  theme(legend.position = "none")

plot_grid(
  plot_grid(p_TPR, p_FDR, 
            ncol = 1),
  legend, 
  ncol = 1, rel_heights = c(2, 0.1))

ggsave("figures/figure_5.png", width = 15, height = 8, dpi = "retina", units = "cm")
