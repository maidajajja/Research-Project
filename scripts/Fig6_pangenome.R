suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Pan-genome composition data from Panaroo results
pg <- data.frame(
  Category = factor(c("Core
(>=99%)", "Soft Core
(95-99%)", "Shell
(15-95%)", "Cloud
(<15%)"),
                    levels=c("Core
(>=99%)", "Soft Core
(95-99%)", "Shell
(15-95%)", "Cloud
(<15%)")),
  Genes = c(3154, 918, 1604, 18023),
  Fill = c("#1A5276","#2E86C1","#E67E22","#C0392B")
)
pg$Pct <- round(pg$Genes / sum(pg$Genes) * 100, 1)
pg$Label <- paste0(formatC(pg$Genes, format="d", big.mark=","), "
(", pg$Pct, "%)")

# Panel A: bar chart of gene categories
fig6a <- ggplot(pg, aes(x=Category, y=Genes, fill=Fill)) +
  geom_col(width=0.6, colour="white", linewidth=0.3) +
  geom_text(aes(label=Label), vjust=-0.3, size=3.4, fontface="bold", colour="grey20", lineheight=0.95) +
  scale_fill_identity() +
  scale_y_continuous(expand=expansion(mult=c(0,0.18)), labels=scales::comma) +
  labs(title=expression(italic("K. pneumoniae")~"Pan-genome Composition"),
       subtitle="Total pan-genome: 23,699 genes (Panaroo v1.3.4, n = 234)",
       x="Gene Category", y="Number of Genes") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,15,10,10))

# Panel B: gene frequency histogram from Rtab
rtab_path <- "/scratch/users/k22017808/KP_Research_Project/06_Pangenome/gene_presence_absence.Rtab"
rtab <- read.table(rtab_path, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
gene_freq <- rowSums(rtab[,-1])
freq_df <- data.frame(freq=gene_freq / 234 * 100)

fig6b <- ggplot(freq_df, aes(x=freq)) +
  geom_histogram(binwidth=5, fill="#2E86C1", colour="white", linewidth=0.2) +
  scale_x_continuous(breaks=seq(0,100,10), labels=function(x) paste0(x,"%")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08)), labels=scales::comma) +
  labs(title="Pan-genome Gene Frequency Distribution",
       subtitle="Frequency of each gene across 234 isolates",
       x="Gene Frequency (% isolates)", y="Number of Genes") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig6a_pangenome_composition.png"), fig6a, width=8, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6b_pangenome_frequency.png"), fig6b, width=8, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6a_pangenome_composition.pdf"), fig6a, width=8, height=6)
ggsave(file.path(OUT,"Fig6b_pangenome_frequency.pdf"), fig6b, width=8, height=6)
message("Fig6 saved")
