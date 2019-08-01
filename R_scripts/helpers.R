facing_trees <- function(tree_cor, tree_phy, data, col_cor, col_phy, p_thresh = 0.01) {
  
  ## Keep only taxa found by one or the other method
  data <- filter(data, {{col_cor}} <= p_thresh | {{col_phy}} <= p_thresh) %>% 
    mutate(Status = case_when(
      {{col_cor}} <= p_thresh & {{col_phy}} <= p_thresh ~ "Both",
      {{col_cor}} <= p_thresh                           ~ "Correlation",
      {{col_phy}} <= p_thresh                           ~ "Phylogeny",
    ))
  
  my_palette <- setNames(c("#C77CFF", "#F8766D", "#349E34"), c("Correlation", "Phylogeny", "Both"))
  
  p1 <- 
    tree_cor %>%
    ggtree(branch.length="none", color = "grey30") %<+%
    data +
    geom_tippoint(aes(size = -log10({{col_cor}}), color = Status), alpha = 0.5) +
    geom_tiplab(aes(color = Status), vjust = 0, hjust = 1, 
                size = 2.5, fontface = "bold") +
    # scale_color_viridis_d(direction = -1) +
    scale_color_manual(values = my_palette) + 
    theme(legend.position = "bottom")
  
  p2 <-
    tree_phy %>% 
    ggtree(color = "grey30") %<+%
    data +
    geom_tippoint(aes(size = -log10({{col_phy}}), color = Status), alpha = 0.5) +
    geom_tiplab(aes(color = Status), vjust = 0, hjust = 0,
                size = 2.5, fontface = "bold") +
    # scale_color_viridis_d(direction = -1) +
    scale_color_manual(values = my_palette) + 
    theme(legend.position = "bottom") + 
    scale_x_reverse()

  plot_grid(p1, p2, 
            ncol = 2, 
            # rel_heights = c(1, 0.1, 1),
            labels = c("A", "B"), 
            label_x = c(0, 0.9), 
            label_y = c(1, 1))
  
}

my_plot <- function(otus, data, env) {
  data %>% rownames_to_column(var = "OTU") %>% gather(-OTU, key = "Sample", value = "Abundance") %>% 
    mutate(food_type = str_sub(Sample, 1, 2)) %>% filter(OTU %in% otus) %>% 
    ggplot(aes(x = {{env}}, y = Abundance, fill = {{env}})) + 
    geom_boxplot(alpha = 0.5) + 
    geom_point(aes(color = {{env}})) + 
    facet_wrap(~OTU, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 90), 
          legend.position = "bottom")
}
