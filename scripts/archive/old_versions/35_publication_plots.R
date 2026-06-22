library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ============================================================
# READ DATA
# ============================================================
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

# Get key columns
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]

kleb$ST <- gsub("-.*LV", "", kleb[[st_col]])
kleb$virulence_score <- as.numeric(kleb[[vir_col]])
kleb$resistance_score <- as.numeric(kleb[[res_col]])

# Major ST colours
major_colours <- c(
  "ST23" = "#00BCD4", "ST86" = "#1565C0", "ST11" = "#FF6F00",
  "ST258" = "#E65100", "ST65" = "#6A1B9A", "ST29" = "#AD1457",
  "ST512" = "#2E7D32", "ST76" = "#00838F", "ST60" = "#558B2F",
  "ST380" = "#795548", "Other" = "#BDBDBD"
)

# ============================================================
# FIGURE 2 — ST DISTRIBUTION
# ============================================================
st_counts <- kleb %>%
  group_by(ST) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(fill_colour = ifelse(ST %in% names(major_colours), ST, "Other"))

# Keep top 20 STs by name, group rest as Other
top_sts <- st_counts %>% head(20)
other_count <- sum(st_counts$n[21:nrow(st_counts)])
other_row <- data.frame(ST = "Other STs", n = other_count, fill_colour = "Other")
plot_data <- bind_rows(top_sts, other_row)
plot_data$ST <- factor(plot_data$ST, levels = plot_data$ST)

fig2 <- ggplot(plot_data, aes(x = ST, y = n, fill = fill_colour)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = major_colours) +
  geom_text(aes(label = n), vjust = -0.4, size = 3, fontface = "bold") +
  labs(
    title = expression(paste("Sequence Type Distribution of 234 ", italic("K. pneumoniae"), " Isolates")),
    x = "Sequence Type",
    y = "Number of Isolates (n)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig2_ST_distribution.png",
       plot = fig2, width = 12, height = 6, dpi = 300, bg = "white")
cat("Figure 2 saved\n")

# ============================================================
# FIGURE 3 — VIRULENCE AND RESISTANCE SCORES
# ============================================================

# 3A - Virulence score
vir_counts <- kleb %>%
  group_by(virulence_score) %>%
  summarise(n = n()) %>%
  mutate(pct = round(n/sum(n)*100, 1))

fig3a <- ggplot(vir_counts, aes(x = factor(virulence_score), y = n, 
                                 fill = factor(virulence_score))) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c("0"="#FFFDE7","1"="#FFF176","2"="#FFD600",
                                "3"="#FF6F00","4"="#7B1FA2","5"="#4A148C")) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")), 
            vjust = -0.3, size = 3, fontface = "bold") +
  labs(
    title = "Virulence Score Distribution",
    x = "Virulence Score (0–5)",
    y = "Number of Isolates (n)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  ) +
  ylim(0, max(vir_counts$n) * 1.15)

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig3a_virulence_distribution.png",
       plot = fig3a, width = 7, height = 6, dpi = 300, bg = "white")

# 3B - Resistance score
res_counts <- kleb %>%
  group_by(resistance_score) %>%
  summarise(n = n()) %>%
  mutate(pct = round(n/sum(n)*100, 1))

fig3b <- ggplot(res_counts, aes(x = factor(resistance_score), y = n,
                                 fill = factor(resistance_score))) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c("0"="#1a9850","1"="#fdae61","2"="#f46d43","3"="#d73027")) +
  geom_text(aes(label = paste0(n, "\n(", pct, "%)")),
            vjust = -0.3, size = 3, fontface = "bold") +
  labs(
    title = "Resistance Score Distribution",
    x = "Resistance Score (0–3)",
    y = "Number of Isolates (n)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  ) +
  ylim(0, max(res_counts$n) * 1.15)

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig3b_resistance_distribution.png",
       plot = fig3b, width = 7, height = 6, dpi = 300, bg = "white")
cat("Figure 3 saved\n")

# ============================================================
# FIGURE 4 — VIRULENCE VS RESISTANCE SCATTER
# ============================================================
kleb$convergent <- ifelse(kleb$virulence_score >= 4 & kleb$resistance_score >= 2,
                          "Convergent", "Non-convergent")
kleb$ST_label <- ifelse(kleb$convergent == "Convergent", kleb$ST, "")

fig4 <- ggplot(kleb, aes(x = resistance_score, y = virulence_score, 
                          colour = convergent, shape = convergent)) +
  geom_jitter(size = 3, alpha = 0.8, width = 0.12, height = 0.12, stroke = 0.5) +
  scale_colour_manual(values = c("Convergent" = "#d73027", "Non-convergent" = "#4575b4"),
                      name = "") +
  scale_shape_manual(values = c("Convergent" = 17, "Non-convergent" = 16), name = "") +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = 3.5, ymax = 5.5,
           fill = NA, colour = "#d73027", linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = 2.5, y = 5.6, label = "Convergent zone\n(n=20)",
           colour = "#d73027", size = 3.5, fontface = "bold") +
  scale_x_continuous(breaks = 0:3) +
  scale_y_continuous(breaks = 0:5) +
  labs(
    title = expression(paste("Virulence vs Resistance Scores in 234 ", italic("K. pneumoniae"), " Isolates")),
    x = "Resistance Score (0–3)",
    y = "Virulence Score (0–5)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3)
  )

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig4_virulence_vs_resistance.png",
       plot = fig4, width = 8, height = 8, dpi = 300, bg = "white")
cat("Figure 4 saved\n")

# ============================================================
# FIGURE 5 — PAN-GENOME COMPOSITION
# ============================================================
categories <- c("Core\n(≥99% isolates)", "Soft Core\n(95–99%)", 
                "Shell\n(15–95%)", "Cloud\n(<15%)")
counts <- c(3154, 918, 1604, 18023)
total <- sum(counts)
pcts <- round(counts/total*100, 1)

df_pan <- data.frame(
  category = factor(categories, levels = rev(categories)),
  count = counts,
  pct = pcts
)

pan_colours <- c("#2166ac", "#74add1", "#fdae61", "#d73027")
names(pan_colours) <- categories

fig5 <- ggplot(df_pan, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.4, width = 0.6) +
  scale_fill_manual(values = c(
    "Core\n(≥99% isolates)" = "#2166ac",
    "Soft Core\n(95–99%)" = "#74add1",
    "Shell\n(15–95%)" = "#fdae61",
    "Cloud\n(<15%)" = "#d73027"
  )) +
  geom_text(aes(label = paste0(format(count, big.mark=","), "\n(", pct, "%)")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(
    title = expression(paste("Pan-genome Composition of 234 ", italic("K. pneumoniae"), " Isolates")),
    subtitle = paste0("Total pan-genome: ", format(total, big.mark=","), " genes"),
    x = "Gene Category",
    y = "Number of Genes"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  ) +
  ylim(0, max(counts) * 1.12)

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig5_pangenome.png",
       plot = fig5, width = 10, height = 7, dpi = 300, bg = "white")
cat("Figure 5 saved\n")

# ============================================================
# FIGURE 6 — TOP AMR GENES PREVALENCE
# ============================================================
amr <- read.table("/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder/amrfinder_all.tsv",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", fill=TRUE)

gene_col <- "Element.symbol"
if(gene_col %in% colnames(amr)) {
  gene_counts <- amr %>%
    group_by(.data[[gene_col]]) %>%
    summarise(n_isolates = n_distinct(Name)) %>%
    arrange(desc(n_isolates)) %>%
    head(20) %>%
    mutate(prevalence = round(n_isolates/234*100, 1))
  
  gene_counts[[gene_col]] <- factor(gene_counts[[gene_col]], 
                                     levels = rev(gene_counts[[gene_col]]))
  
  fig6 <- ggplot(gene_counts, aes(x = .data[[gene_col]], y = prevalence)) +
    geom_bar(stat = "identity", fill = "#d73027", colour = "black", 
             linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = paste0(n_isolates, "\n(", prevalence, "%)")),
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
    xlim(0, 110)
  
  ggsave("/scratch/users/k22017808/KP_Research_Project/plots/Fig6_AMR_prevalence.png",
         plot = fig6, width = 10, height = 8, dpi = 300, bg = "white")
  cat("Figure 6 saved\n")
}

cat("\nAll publication figures complete!\n")
