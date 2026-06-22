library(ggplot2)

# Pan-genome composition
categories <- c("Core\n(≥99%)", "Soft Core\n(95-99%)", "Shell\n(15-95%)", "Cloud\n(<15%)")
counts <- c(3154, 918, 1604, 18023)
total <- sum(counts)
percentages <- round(counts/total*100, 1)

df <- data.frame(
  category = factor(categories, levels = categories),
  count = counts,
  percentage = percentages
)

colours <- c("#2166ac", "#74add1", "#fdae61", "#d73027")

p1 <- ggplot(df, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, colour = "white", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colours) +
  labs(title = "Pan-genome Composition of 234 K. pneumoniae Isolates",
       subtitle = paste0("Total pan-genome: ", total, " genes"),
       fill = "Gene Category") +
  theme_void() +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right") +
  geom_text(aes(label = paste0(percentage, "%\n(", count, ")")),
            position = position_stack(vjust = 0.5), size = 3.5, colour = "white", fontface = "bold")

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/pangenome_composition.png",
       plot = p1, width = 10, height = 7, dpi = 300)

# Bar chart version
p2 <- ggplot(df, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = colours) +
  geom_text(aes(label = paste0(count, "\n(", percentage, "%)")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(title = "Pan-genome Composition of 234 K. pneumoniae Isolates",
       x = "Gene Category", y = "Number of Genes") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "none") +
  ylim(0, max(counts) * 1.1)

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/pangenome_barchart.png",
       plot = p2, width = 10, height = 7, dpi = 300)

cat("Pan-genome plots saved\n")
