mytheme <-
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text = element_text(size = 8), 
                 axis.title = element_text(size = 12),
                 legend.position = "none")


color_values <- c(
  "Correlation" = "#C77CFF",
  "Bootstrap" = "#00BFC4",
  "Random Correlation" = "#7CAE00",
  "Phylogeny" = "#F8766D",
  "Taxonomy" = "#F8766D",
  "Random Phylogeny" = "#FFA500",
  "Random Taxonomy" = "#FFA500",
  "BH" = "#4169E1"
)

size_values <- c(
  "Correlation" = 4,
  "Bootstrap" = 1,
  "Random Correlation" = 1,
  "Phylogeny" = 4,
  "Taxonomy" = 4,
  "Random Phylogeny" = 1,
  "Random Taxonomy" = 1
)

alpha_values <- c(
  "Correlation" = 0.8,
  "Bootstrap" = .4,
  "Random Correlation" = .4,
  "Phylogeny" = 0.8,
  "Taxonomy" = 0.8,
  "Random Phylogeny" = .4,
  "Random Taxonomy" = .4
)

shape_values <- c(
  "Correlation" = 17,
  "Bootstrap" = 2,
  "Random Correlation" = 6,
  "Phylogeny" = 16,
  "Taxonomy" = 16,
  "Random Phylogeny" = 1,
  "Random Taxonomy" = 1
)

linetype_values <- c(
  "Correlation" = "solid",
  "Taxonomy" = "dotdash",
  "Random Correlation" = "dashed",
  "Random Taxonomy" = "twodash",
  "BH" = "dotted"
)