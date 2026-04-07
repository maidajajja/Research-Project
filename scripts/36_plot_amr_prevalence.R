library(ggplot2)
library(dplyr)

amr <- read.table("/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder/amrfinder_all.tsv",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", fill=TRUE,
                  check.names=TRUE)

# Print actual column names
cat("Columns:", paste(colnames(amr)[1:10], collapse=", "), "\n")

# Use check.names=TRUE so spaces become dots
gene_col <- "Element.symbol"
name_col <- colnames(amr)[1]  # first column is sample name
cat("Name column:", name_col, "\n")

gene_counts <- amr %>%
  group_by(.data[[gene_col]]) %>%
  summarise(n_isolates = n_distinct(.data[[name_col]])) %>%
  arrange(desc(n_isolates)) %>%
  head(20) %>%
  mutate(prevalence = round(n_isolates/234*100, 1))

gene_counts[[gene_col]] <- factor(gene_counts[[gene_col]],
                                   levels = rev(gene_counts[[gene_col]]))

fig6 <- ggplot(gene_counts, aes(x = .data[[gene_col]], y = prevalence)) +
  geom_bar(stat = "identity", fill = "#d73027", colour = "black",
           linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = paste0(n_isolates, " (", prevalence, "%)")),
            hjust = -0.1, size = 3, fontface = "bold") +
  coord_flip() +
  labs(
    title = expression(paste("Prevalence of AMR Genes in 234 ", italic("K. pneumoniae"), " Isolates")),
    x = "AMR Gene",
    y = "Prevalence (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
  ) +
  scale_y_continuous(limits = c(0, 115))

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig6_AMR_prevalence.png",
       plot = fig6, width = 10, height = 8, dpi = 300, bg = "white")
cat("Figure 6 saved\n")
