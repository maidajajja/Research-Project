library(ggplot2)
library(dplyr)

df <- read.csv("ST23.results.csv", stringsAsFactors = FALSE)

finite_max <- max(df$Odds_ratio[is.finite(df$Odds_ratio)], na.rm = TRUE)
df$Odds_ratio_plot <- ifelse(is.infinite(df$Odds_ratio), finite_max * 1.2, df$Odds_ratio)
df$log_p <- -log10(pmax(df$Bonferroni_p, 1e-300))

named_genes <- c("glxR", "ybbY", "allA", "allB", "allC", "allD", "allH", "allR", "allS")
df$is_named <- df$Gene %in% named_genes
df$significant <- df$Bonferroni_p < 0.05

# Get the shared coordinate of the named-gene point (they all sit at the same spot)
named_point <- df %>% filter(is_named) %>% slice(1)
label_text <- paste(sort(named_genes), collapse = ", ")

p <- ggplot(df, aes(x = Odds_ratio_plot, y = log_p)) +
  geom_point(aes(color = significant), alpha = 0.4, size = 1.5, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "#4477AA", "FALSE" = "grey70")) +
  geom_point(data = subset(df, is_named), color = "#CC3311", size = 3) +
  annotate("segment",
           x = named_point$Odds_ratio_plot, y = named_point$log_p,
           xend = named_point$Odds_ratio_plot * 0.55, yend = named_point$log_p + 6,
           color = "grey30", linewidth = 0.4) +
  annotate("label",
           x = named_point$Odds_ratio_plot * 0.55, y = named_point$log_p + 6.5,
           label = label_text, hjust = 1, vjust = 0, size = 3.3, fontface = "italic",
           fill = "white", label.size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "Odds ratio (capped for infinite values)",
       y = expression(-log[10]~"(Bonferroni-corrected p-value)"),
       title = "Accessory genes associated with ST23 membership") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 13, face = "bold")) +
  coord_cartesian(clip = "off")

ggsave("FigS5_scoary_volcano_v2.png", p, width = 9, height = 6.5, dpi = 600)
ggsave("FigS5_scoary_volcano_v2.pdf", p, width = 9, height = 6.5)

cat("Done.\n")
