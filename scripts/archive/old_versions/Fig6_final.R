suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# ── Fig6a: Pan-genome composition bar chart ──────────────────────────────────
pg <- data.frame(
  Category = factor(c("Core\n(\u226599%)", "Soft core\n(95\u201399%)",
                       "Shell\n(15\u201395%)", "Cloud\n(<15%)"),
                    levels=c("Core\n(\u226599%)", "Soft core\n(95\u201399%)",
                             "Shell\n(15\u201395%)", "Cloud\n(<15%)")),
  Genes = c(3154, 918, 1604, 18023)
)
pg$Pct   <- round(pg$Genes / sum(pg$Genes) * 100, 1)
pg$Label <- paste0(formatC(pg$Genes, format="d", big.mark=","),
                   "\n(", pg$Pct, "%)")

# Wong CBF palette mapped to pan-genome categories
cbf_pg <- c("#0072B2","#56B4E9","#E69F00","#D55E00")

fig6a <- ggplot(pg, aes(x=Category, y=Genes, fill=Category)) +
  geom_col(width=0.6, colour="white", linewidth=0.3) +
  geom_text(aes(label=Label), vjust=-0.3, size=3.5,
            fontface="bold", colour="grey20", lineheight=0.95) +
  scale_fill_manual(values=setNames(cbf_pg, levels(pg$Category)),
                    guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18)),
                     labels=scales::comma) +
  labs(x="Gene category", y="Number of genes") +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig6a_pangenome_composition.png"), fig6a,
       width=7, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6a_pangenome_composition.pdf"), fig6a,
       width=7, height=6)

# ── Fig6b: Gene frequency distribution ──────────────────────────────────────
rtab_path <- "/scratch/users/k22017808/KP_Research_Project/06_Pangenome/gene_presence_absence.Rtab"
rtab <- read.table(rtab_path, sep="\t", header=TRUE,
                   stringsAsFactors=FALSE, check.names=FALSE)
gene_freq <- rowSums(rtab[,-1])
n_isolates <- ncol(rtab) - 1
freq_df <- data.frame(freq = gene_freq / n_isolates * 100)

fig6b <- ggplot(freq_df, aes(x=freq)) +
  geom_histogram(binwidth=5, fill="#0072B2", colour="white",
                 linewidth=0.2, boundary=0) +
  # threshold lines
  geom_vline(xintercept=15, linetype="dashed",
             colour="#E69F00", linewidth=0.8) +
  geom_vline(xintercept=95, linetype="dashed",
             colour="#0072B2", linewidth=0.8) +
  # annotations
  annotate("text", x=7, y=Inf, label="Cloud\n(<15%)",
           vjust=1.4, hjust=0.5, size=3.2, colour="#E69F00", fontface="bold") +
  annotate("text", x=97, y=Inf, label="Core\n(\u226599%)",
           vjust=1.4, hjust=0, size=3.2, colour="#0072B2", fontface="bold") +
  scale_x_continuous(breaks=seq(0,100,10),
                     labels=function(x) paste0(x,"%")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08)),
                     labels=scales::comma) +
  labs(x="Gene frequency (% isolates)", y="Number of genes") +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig6b_pangenome_frequency.png"), fig6b,
       width=8, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6b_pangenome_frequency.pdf"), fig6b,
       width=8, height=6)

message("Fig6 final saved")
