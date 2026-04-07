library(ggplot2)
library(dplyr)
library(reshape2)

# Read AMRFinder results
amr <- read.table("/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder/amrfinder_all.tsv",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", fill=TRUE)

# Get top 20 most common genes
gene_col <- "Element.symbol"
sample_col <- "Name"

if(gene_col %in% colnames(amr) & sample_col %in% colnames(amr)) {
  gene_counts <- amr %>%
    group_by(.data[[gene_col]]) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    head(20)
  
  top_genes <- gene_counts[[gene_col]]
  
  amr_sub <- amr %>% filter(.data[[gene_col]] %in% top_genes)
  
  # Create presence/absence matrix
  amr_matrix <- amr_sub %>%
    mutate(present = 1) %>%
    group_by(.data[[sample_col]], .data[[gene_col]]) %>%
    summarise(present = max(present), .groups = "drop")
  
  p <- ggplot(amr_matrix, aes(x = .data[[gene_col]], 
                               y = .data[[sample_col]], 
                               fill = factor(present))) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_fill_manual(values = c("1" = "#d73027"), na.value = "#f7f7f7") +
    labs(title = "AMR Gene Distribution Across 234 K. pneumoniae Isolates",
         subtitle = "Top 20 most prevalent resistance genes",
         x = "AMR Gene", y = "Isolate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "none")
  
  ggsave("/scratch/users/k22017808/KP_Research_Project/plots/amr_heatmap.png",
         plot = p, width = 12, height = 10, dpi = 300)
  cat("AMR heatmap saved\n")
} else {
  cat("Column names:", paste(colnames(amr)[1:10], collapse=", "), "\n")
}
