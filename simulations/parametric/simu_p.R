library(curatedMetagenomicData)
library(correlationtree)
library(StructFDR)
library(tidyverse)
library(phyloseq)
library(janitor)
library(cowplot)
library(evabic)
library(furrr)
library(yatah)
library(glue)

options("future.fork.enable" = TRUE)
plan(multiprocess)


#### Data ####

N_sample <- 100
samples <- paste0("S", c(paste0("0", c(paste0("0", 1:9), 10:99)), 100))

Notus <- 400

# saveRDS(Notus, "simulations/parametric/simus_p-notus.rds")
# Notus <- readRDS("simulations/parametric/simus_p-notus.rds")

nb_mean <- 10000
nb_size <- 25

data("alcohol")
data <- alcohol$X

## Abundance

otus <- 
  data %>% 
  rowSums() %>% 
  sort(decreasing = TRUE) %>% 
  head(n = Notus) %>% 
  names()

data_filtered <- data[otus, ]

otu_table_filtered <- alcohol$otu.name[otus, ]

## Fit

fit <- dirmult(t(data_filtered), epsilon = 10^(-5), trace = FALSE)
gamma <- fit$gamma
names(gamma) <- otus

## Simulation

param_nb <- rnbinom(N_sample, mu = nb_mean, size = nb_size)
param_dirich <- rerun(N_sample, rdirichlet(1, alpha = gamma))

data_sim <-
  map2_dfc(param_nb, param_dirich, rmultinom, n = 1) %>% 
  as.matrix() %>% 
  `colnames<-`(samples) %>% 
  `rownames<-`(otus)

## Trees

tree_cor <- correlation_tree(data_sim, matrix = TRUE, method = "spearman")


#### Phylogenetic tree #### 

tree_phy <- prune_taxa(otus, alcohol$tree)
tree_phy <- multi2di(tree_phy)
tree_phy$node.label <- NULL
tree_phy$edge.length <- 
  mean_lineage_length(tree_cor) * tree_phy$edge.length / 
  mean_lineage_length(tree_phy)



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

B <- 20 # In the paper, B = 100

my_TreeFDR <- partial(TreeFDR2, B = B, q.cutoff = 0.5,
                      test.func = test.func.wt, perm.func = perm.func)


## Parameters

fc <- c(5, 10, 15, 20)
nH1 <- c(2, 10, 25, 40) # In the paper, nH1 = c(2, 5, 10, 15, 25, 40)
repl <- 4 # In the paper, repl > 600 

set.seed(42)

## Initialisation

time0 <- format(Sys.time(), "%y-%m-%d_%H-%M")

df_simus <-
  crossing(fc, nH1) %>% 
  rerun(repl, .) %>% 
  bind_rows() %>% 
  mutate(time = time0, B = B) %>% 
  rowid_to_column("ID") %>% 
  mutate(ind_samp = rerun(n(), sample(N_sample, round(N_sample / 2))),
         samp_lgl = map(ind_samp, ~ seq(N_sample) %in% .),
         otus_diffs = map(nH1, sample, x = otus),
         ind_otus = map(otus_diffs, ~ which(otus %in% .)),
         newdata = pmap(list(fc = fc, rows = ind_otus, cols = ind_samp),
                        apply_fc, X = data_sim), 
         tree_cor = future_map(newdata, correlation_tree, matrix = TRUE, method = "spearman"),
         tree_randphy = rerun(n(), shuffle_tiplabels(tree_phy)),
         tree_randcor = future_map(tree_cor, shuffle_tiplabels))


## Analysis

df_simus <- 
  df_simus %>% # Could take several hours
  mutate(fdrobj_cor = future_pmap(list(X = newdata, Y = samp_lgl, tree = tree_cor), my_TreeFDR),
         fdrobj_randcor = future_pmap(list(X = newdata, Y = samp_lgl, tree = tree_randcor), my_TreeFDR),
         fdrobj_phy = future_pmap(list(X = newdata, Y = samp_lgl), my_TreeFDR, tree = tree_phy),
         fdrobj_randphy = future_pmap(list(X = newdata, Y = samp_lgl, tree = tree_randphy), my_TreeFDR)) %>% 
  select(-newdata, -ind_samp, -ind_otus, -samp_lgl, -tree_cor, -tree_randphy, -tree_randcor)

## Too big to be committed
# saveRDS(df_simus, "simulations/parametric/simus_p-df_simus.rds")
# df_simus <- readRDS("simulations/parametric/simus_p-df_simus.rds")

df_gathered <-
  df_simus %>% 
  gather(method, fdr_obj, -ID, -fc, -nH1, -time, -B, -otus_diffs) %>% 
  mutate(method = str_remove_all(method, "fdrobj_"))

## Too big to be committed
# saveRDS(df_gathered, "simulations/parametric/simus_p-df_gathered.rds")
# df_gathered <- readRDS("simulations/parametric/simus_p-df_gathered.rds")

## To save the results in several <100 MB files (that can be commited)
# allgroups <- letters[1:12]
# 
# groups <- sort(rep_along(seq_len(nrow(df_gathered)), allgroups))
# allgroups %>%
#   map(~ which(groups == .)) %>%
#   map(~ slice(df_gathered, .)) %>%
#   walk2(allgroups,
#         ~ saveRDS(.x, paste0("simulations/parametric/simus_p-df_gathered_part",
#                              .y, ".rds")))
# df_gathered <-
#   paste0("simulations/parametric/simus_p-df_gathered_part", allgroups, ".rds") %>%
#   map(readRDS) %>%
#   reduce(rbind)

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

source("figures/theme.R")

labels <- c("5" = "fc = 5", "10" = "fc = 10", "15" = "fc = 15", "20" = "fc = 20")

## Smoothing

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

ggsave("simulations/parametric/simu_p-smoothing.png", width = 15, height = 5, dpi = "retina")

## TPR and FDR 

df_ebc <- 
  df_eval %>% 
  arrange(nH1) %>% 
  mutate(nH1 = as_factor(nH1)) %>% 
  mutate(method = factor(method, levels = c("bh", "cor", "phy", "randcor", "randphy"), 
                         labels = c("BH", "Correlation", "Phylogeny",
                                    "Random Correlation", "Random Phylogeny")))

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
  scale_linetype_manual(values = linetype_values) +
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
  scale_linetype_manual(values = linetype_values) +
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

ggsave("simulations/parametric/simu_p-ebc.png", width = 15, height = 8, dpi = "retina")
