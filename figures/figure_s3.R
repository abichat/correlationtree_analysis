library(tidyverse)
library(evabic)

source("figures/theme.R")

Notus <- readRDS("simulations/parametric/simus_p-notus.rds")
df_gathered <-
  list.files("simulations/parametric/", 
             pattern = "simus_p-df_gathered_part[a-z].rds") %>%
  paste0("simulations/parametric/", .) %>% 
  map(readRDS) %>%
  reduce(rbind)

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
  select(fc, nH1, otus_diffs, method, B, detected, k, rho, smoothing_mean)

df_treefdr <- 
  df_gathered %>% 
  drop_na(fdr_obj) %>% 
  mutate(detected = map(fdr_obj, ~ names(.$p.unadj)[.$p.adj < 0.05])) %>% 
  mutate(k = map_dbl(fdr_obj, "k"), 
         rho = map_dbl(fdr_obj, "rho")) %>% 
  mutate(smoothing_mean = map_dbl(fdr_obj, ~ mean(abs(.$z.adj - .$z.unadj)))) %>% 
  select(fc, nH1, otus_diffs, method, B, detected, k, rho, smoothing_mean)

df_eval <-
  rbind(df_treefdr, df_bh) %>% 
  mutate(pi0 = 100 * (Notus - nH1) / Notus, 
         tidyebc = map2(detected, otus_diffs, ebc_tidy, m = Notus,
                        measures = c("BACC", "ACC", "TPR", "FDR", "F1"))) %>% 
  unnest(tidyebc) %>% 
  mutate(BACC = ifelse(is.nan(BACC), 0, BACC)) %>% 
  mutate(FDR = ifelse(is.nan(FDR), 0, FDR)) %>% 
  select(-otus_diffs, -detected)

df_smoothing <- 
  df_eval %>% 
  arrange(nH1) %>% 
  mutate(nH1 = as_factor(nH1)) %>% 
  mutate(method = factor(method, levels = c("bh", "cor", "phy", "randcor", "randphy"), 
                         labels = c("BH", "Correlation", "Phylogeny",
                                    "Random Correlation", "Random Phylogeny")),
         method = fct_rev(method)) %>% 
  filter(method != "bh")

ggplot(df_smoothing) +
  aes(x = smoothing_mean, fill = method, color = method) +
  geom_density(alpha = 0.7, size = 1, adjust = 1) +
  scale_x_log10(breaks = 10^(-5*0:5)) +
  scale_color_manual(values = color_values, name = "Method",
                     aesthetics = c("color", "fill"), breaks = rev) + 
  labs(x = "Mean z-smoothing", y = "Density") +
  theme_minimal() +
  theme(legend.position = c(0.15, 0.8), 
        legend.justification = c(0, 1), 
        legend.background = element_blank(), 
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 23))

ggsave("figures/figure_s3.png", width = 15, height = 5, dpi = "retina")

