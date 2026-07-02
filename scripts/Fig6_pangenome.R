suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

pg <- data.frame(
  Category = factor(c("Core (>=99%)", "Soft core (95-99%)", "Shell (15-95%)", "Cloud (<15%)"),
                    levels=c("Core (>=99%)", "Soft core (95-99%)", "Shell (15-95%)", "Cloud (<15%)")),
  Genes = c(3769, 368, 1558, 8356),
  Fill = c("#1A5276","#2E86C1","#E67E22","#C0392B")
)
pg$Pct <- round(pg$Genes / sum(pg$Genes) * 100, 1)
pg$Label <- paste0(formatC(pg$Genes, format="d", big.mark=","), "\n(", pg$Pct, "%)")

fig6a <- ggplot(pg, aes(x=Category, y=Genes, fill=Fill)) +
  geom_col(width=0.55, colour="white", linewidth=0.3) +
  geom_text(aes(label=Label), vjust=-0.3, size=4.5, fontface="bold", colour="grey20", lineheight=0.95) +
  scale_fill_identity() +
  scale_y_continuous(expand=expansion(mult=c(0,0.18)), labels=scales::comma) +
  labs(
       x="Gene Category", y="Number of Genes") +
  theme_classic(base_size=15) +
  theme(plot.subtitle=element_text(size=12, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=14),
        axis.text.x=element_text(size=13, colour="grey20"),
        axis.text.y=element_text(size=13, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,20,10,10))

rtab_path <- "/scratch/users/k22017808/KP_Research_Project/06_Pangenome_final229/gene_presence_absence.Rtab"
rtab <- read.table(rtab_path, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
gene_freq <- rowSums(rtab[,-1])
freq_df <- data.frame(freq=gene_freq / 229 * 100)

fig6b <- ggplot(freq_df, aes(x=freq)) +
  geom_histogram(binwidth=5, fill="#2E86C1", colour="white", linewidth=0.2) +
  geom_vline(xintercept=15, linetype="dashed", colour="#E67E22", linewidth=0.8) +
  geom_vline(xintercept=99, linetype="dashed", colour="#1A5276", linewidth=0.8) +
  annotate("text", x=16, y=14000, label="Cloud\n(<15%)", colour="#E67E22", hjust=0, size=4.5, fontface="bold") +
  annotate("text", x=85, y=14000, label="Core\n(>=99%)", colour="#1A5276", hjust=0, size=4.5, fontface="bold") +
  scale_x_continuous(breaks=seq(0,100,10), labels=function(x) paste0(x,"%")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08)), labels=scales::comma) +
  labs(
       x="Gene Frequency (% isolates)", y="Number of Genes") +
  theme_classic(base_size=15) +
  theme(plot.subtitle=element_text(size=12, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=14),
        axis.text=element_text(size=13, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,20,10,10))

ggsave(file.path(OUT,"Fig6a_pangenome_composition.png"), fig6a, width=9, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6b_pangenome_frequency.png"), fig6b, width=9, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6a_pangenome_composition.pdf"), fig6a, width=9, height=7)
ggsave(file.path(OUT,"Fig6b_pangenome_frequency.pdf"), fig6b, width=9, height=7)
message("Fig6 saved")
