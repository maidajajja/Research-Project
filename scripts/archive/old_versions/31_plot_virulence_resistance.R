library(ggplot2)
library(dplyr)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

vir_col <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]

kleb$virulence_score <- as.numeric(kleb[[vir_col]])
kleb$resistance_score <- as.numeric(kleb[[res_col]])
kleb$ST <- gsub("-.*LV", "", kleb[[st_col]])

# Plot 1 - Virulence score distribution
p1 <- ggplot(kleb, aes(x = factor(virulence_score), fill = factor(virulence_score))) +
  geom_bar(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("0"="#FFFDE7","1"="#FFF176","2"="#FFD600",
                                "3"="#FF6F00","4"="#7B1FA2","5"="#4A148C")) +
  labs(title = "Virulence Score Distribution (Kleborate)",
       subtitle = "234 K. pneumoniae isolates from liver disease patients",
       x = "Virulence Score (0-5)", y = "Number of Isolates", fill = "Score") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "none")

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/virulence_distribution.png",
       plot = p1, width = 8, height = 6, dpi = 300)

# Plot 2 - Resistance score distribution
p2 <- ggplot(kleb, aes(x = factor(resistance_score), fill = factor(resistance_score))) +
  geom_bar(colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("0"="#1a9850","1"="#fdae61","2"="#f46d43","3"="#d73027")) +
  labs(title = "Resistance Score Distribution (Kleborate)",
       subtitle = "234 K. pneumoniae isolates from liver disease patients",
       x = "Resistance Score (0-3)", y = "Number of Isolates", fill = "Score") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "none")

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/resistance_distribution.png",
       plot = p2, width = 8, height = 6, dpi = 300)

# Plot 3 - Virulence vs resistance scatter
kleb$convergent <- ifelse(kleb$virulence_score >= 4 & kleb$resistance_score >= 2,
                          "Convergent (n=20)", "Non-convergent")

p3 <- ggplot(kleb, aes(x = resistance_score, y = virulence_score, colour = convergent)) +
  geom_jitter(size = 3, alpha = 0.7, width = 0.15, height = 0.15) +
  scale_colour_manual(values = c("Convergent (n=20)" = "#d73027",
                                  "Non-convergent" = "#4575b4")) +
  labs(title = "Virulence vs Resistance Score per Isolate",
       subtitle = "Convergent strains carry both high virulence and resistance",
       x = "Resistance Score (0-3)", y = "Virulence Score (0-5)", colour = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom")

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/virulence_vs_resistance.png",
       plot = p3, width = 8, height = 8, dpi = 300)

cat("Virulence and resistance plots saved\n")
