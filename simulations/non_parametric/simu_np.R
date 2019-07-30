library(curatedMetagenomicData)
library(correlationtree)
library(tidyverse)
library(janitor)
library(cowplot)
library(evabic)
library(furrr)
library(yatah)
library(glue)

options("future.fork.enable" = TRUE)
plan(multiprocess)


#### Data ####

exprSet <-
  "BritoIL_2016.metaphlan_bugs_list.stool" %>%
  curatedMetagenomicData(dryrun = FALSE, counts = TRUE) %>%
  mergeData()

## Samples

df_sample <-
  exprSet %>%
  pData() %>%
  as_tibble(rownames = "Sample") %>%
  select(Sample, age_category) %>%
  filter(age_category == "adult") %>%
  mutate(Status = sample(c("A", "B"), size = n(), replace = TRUE),
         age_category = NULL)

Nsamples <- nrow(df_sample)

## Abundance 

df_abund <-
  exprSet %>%
  exprs() %>%
  as_tibble(rownames = "taxonomy") %>%
  filter(is_rank(taxonomy, "genus")) %>%
  mutate(clade = last_clade(taxonomy)) %>%
  select(clade, taxonomy, one_of(df_sample$Sample))

clades_to_keep <-
  df_abund %>%
  select(-taxonomy) %>%
  gather(key = Sample, value = Count, -clade) %>%
  group_by(clade) %>%
  summarise(P = sum(Count > 0), M = mean(Count > 0)) %>%
  arrange(P) %>%
  filter(M > 0.2) %>%
  pull(clade)

df_abund <- filter(df_abund, clade %in% clades_to_keep)

## Trees

tree_cor <-
  df_abund %>%
  select(-taxonomy) %>%
  correlation_tree(method = "spearman")

tree_tax <-
  df_abund %>%
  pull(taxonomy) %>%
  taxtable() %>%
  taxtree(lineage_length = mean_lineage_length(tree_cor))

tree_tax$node.label <- NULL

## Taxa

taxa_all <-
  df_abund %>%
  select(-taxonomy) %>%
  gather(key = Sample, value = Count, -clade) %>% 
  group_by(clade) %>% 
  summarise(m = length(Count[Count > 0])) %>% 
  arrange(desc(m)) %>%
  pull(clade) 

Ntaxa <- length(taxa_all)

taxa_top30 <- taxa_all[1:30]


#### Simulations ####

## Functions

apply_fc <- function(X, fc = 1, rows = 0, cols = 0){
  
  if(length(rows) * length(cols) == 0){
    return(X)
  }
  
  if(max(rows) > nrow(X) | max(cols) > ncol(X)){
    stop("Indexes out of range.")
  }
  
  mat <- matrix(1, nrow = nrow(X), ncol = ncol(X))
  mat[rows, cols] <- fc
  
  X * mat
}

perm.func <- function (X, Y, ...) {
  return(list(X = X, Y = sample(Y)))
}

test.func.wt <- function (X, Y) {
  Y <- as.numeric(factor(Y))
  obj <- apply(X, 1, function(x) {
    p.value <- suppressWarnings(wilcox.test(x ~ Y)$p.value)
    e.sign <- sign(mean(x[Y == 2]) - mean(x[Y == 1]))
    c(p.value, e.sign)
  })
  return(list(p.value=obj[1, ], e.sign=obj[2, ]))
}

B <- 20 # In the papare, B = 100

my_TreeFDR <- partial(TreeFDR2, B = B, q.cutoff = 0.5,
                      test.func = test.func.wt, perm.func = perm.func)

## Formatting

X <-
  df_abund %>%
  select(-taxonomy) %>%
  column_to_rownames("clade") %>%
  as.matrix()

taxa <- rownames(X)

## Parameters

fc <- c(5, 10, 15, 20)
nH1 <- c(2, 5, 10, 15)
repl <- 30 # In the paper, repl > 600 

set.seed(42)

## Initialisation

time0 <- format(Sys.time(), "%y-%m-%d_%H-%M")

df <-
  crossing(fc, nH1) %>% 
  rerun(repl, .) %>% 
  bind_rows() %>% 
  mutate(time = time0, B = B) %>% 
  rowid_to_column("ID") %>% 
  mutate(ind_samp = rerun(n(), sample(Nsamples, round(Nsamples / 2))),
         samp_lgl = map(ind_samp, ~ seq(Nsamples) %in% .),
         taxa_diffs = map(nH1, sample, x = taxa_top30),
         ind_taxa = map(taxa_diffs, ~ which(taxa %in% .)),
         newdata = pmap(list(fc = fc, rows = ind_taxa, cols = ind_samp),
                        apply_fc, X = X),
         tree_cor = future_map(newdata, correlation_tree, matrix = TRUE, method = "spearman"),
         tree_randtax = rerun(n(), shuffle_tiplabels(tree_tax)),
         tree_randcor = future_map(tree_cor, shuffle_tiplabels))

## Analysis

df <- 
  df %>% # Could take several hours
  mutate(fdrobj_cor = future_pmap(list(X = newdata, Y = samp_lgl, tree = tree_cor), my_TreeFDR),
         fdrobj_randcor = future_pmap(list(X = newdata, Y = samp_lgl, tree = tree_randcor), my_TreeFDR),
         fdrobj_tax = future_pmap(list(X = newdata, Y = samp_lgl), my_TreeFDR, tree = tree_tax),
         fdrobj_randtax = future_pmap(list(X = newdata, Y = samp_lgl, tree = tree_randtax), my_TreeFDR)) %>% 
  select(-newdata, -ind_samp, -ind_taxa, -samp_lgl, -tree_cor, -tree_randtax, -tree_randcor)

df_gathered <-
  df %>% 
  gather(method, fdr_obj, -ID, -fc, -nH1, -time, -B, -taxa_diffs) %>% 
  mutate(method = str_remove_all(method, "fdrobj_"))

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
  mutate(k = NA, rho = NA, smoothing = NA, 
         smoothing_med = NA, smoothing_mean = NA, smoothing_Q3 = NA, smoothing_max = NA) %>%
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
  select(-taxa_diffs, -detected)


#### Results ####

df_gathered %>% 
  mutate(null = map_lgl(fdr_obj, is.null),
         null = ifelse(null, "failed", "succeed")) %>%
  tabyl(method, null) %>%
  adorn_totals("row") %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting() %>%
  adorn_ns() %>%
  adorn_title("top", col_name = "")

df_eval %>% 
  group_by(method) %>% 
  summarise_at(vars(BACC, ACC, TPR, FDR), mean)


#### Plots ####

color_values <- c("Correlation" = "#C77CFF", "Taxonomy" = "#F8766D", 
                  "Random Correlation" = "#7CAE00", "Random Taxonomy" = "#FFA500",
                  "BH" = "#4169E1")

linetype_values <- c("Correlation" = "solid", "Taxonomy" = "dotdash", 
                     "Random Correlation" = "dashed", "Random Taxonomy" = "twodash",
                     "BH" = "dotted")

labels <- c("5" = "fc = 5", "10" = "fc = 10", "15" = "fc = 15", "20" = "fc = 20")

## Smoothing

df_smoothing <- 
  df_eval %>% 
  arrange(nH1) %>% 
  mutate(nH1 = as_factor(nH1)) %>% 
  mutate(method = factor(method, levels = c("bh", "cor", "tax", "randcor", "randtax"), 
                         labels = c("BH", "Correlation", "Taxonomy",
                                    "Random Correlation", "Random Taxonomy")),
         method = fct_rev(method)) %>% 
  filter(fc %in% c(5, 10, 15, 20),
         nH1 %in% c(2, 5, 10, 15),
         method != "bh")

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

ggsave("simulations/non_parametric/simu_np-smoothing.png", width = 15, height = 5, dpi = "retina")

## TPR and FDR 

df_ebc <- 
  df_eval %>% 
  arrange(nH1) %>% 
  mutate(nH1 = as_factor(nH1)) %>% 
  mutate(method = factor(method, levels = c("bh", "cor", "tax", "randcor", "randtax"), 
                         labels = c("BH", "Correlation", "Taxonomy",
                                    "Random Correlation", "Random Taxonomy"))) %>% 
  filter(fc %in% c(5, 10, 15, 20),
         nH1 %in% c(2, 5, 10, 15))

# TPR
df_TPR <-
  df_ebc %>% 
  group_by(pi0, fc, method) %>% 
  summarise(mean = mean(TPR), sd = sd(TPR), count = n()) %>% 
  arrange(desc(method)) %>% 
  mutate(infbound = mean - sd/sqrt(count), 
         supbound = mean + sd/sqrt(count))

p_TPR <-
  ggplot(df_TPR) +
  aes(x = pi0, y = mean, color = method) +
  geom_errorbar(aes(ymin = infbound, ymax = supbound, width = 1)) +
  geom_line(aes(linetype = method, group = method)) +
  geom_point(show.legend = FALSE) +
  scale_color_manual(values = color_values) +
  facet_wrap(~ fc, ncol = 4, labeller = labeller(fc = labels)) +
  labs(x = NULL, y = "TPR",
       color = "Method", linetype = "Method") +
  theme_bw() +
  theme(legend.position = "bottom")

legend <- get_legend(p_TPR)

p_TPR <- p_TPR + theme(legend.position = "none")


# FDR
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
  facet_wrap(~ fc, ncol = 4, labeller = labeller(fc = labels)) +
  labs(x = "Proportion of null hypothesis", y = "FDR",
       color = "Method", linetype = "Method") +
  theme_bw() +
  theme(legend.position = "none")

# Combined
plot_grid(
  plot_grid(p_TPR, p_FDR, 
            ncol = 1),
  legend, 
  ncol = 1, rel_heights = c(2, 0.1))

ggsave("simulations/non_parametric/simu_np-ebc.png", width = 15, height = 8, dpi = "retina")
