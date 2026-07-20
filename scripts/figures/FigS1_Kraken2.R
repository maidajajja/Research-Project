suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

kraken <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/kraken2_summary.tsv",
                     sep="\t", header=TRUE, stringsAsFactors=FALSE)
colnames(kraken) <- c("Sample","Klebsiella_pct","Kp_pct","Unclassified_pct")
kraken <- kraken[order(kraken$Kp_pct), ]
kraken$Sample <- factor(kraken$Sample, levels=kraken$Sample)
kraken$Pass <- kraken$Kp_pct >= 80

figS1 <- ggplot(kraken, aes(x=Sample, y=Kp_pct, fill=Pass)) +
  geom_col(width=0.8, colour=NA) +
  geom_hline(yintercept=80, linetype="dashed", colour="#C0392B", linewidth=0.7) +
  annotate("text", x=5, y=82, label="80% threshold", colour="#C0392B",
           size=3.5, fontface="bold", hjust=0) +
  scale_fill_manual(values=c("TRUE"="#2471A3","FALSE"="#E74C3C"),
                    labels=c("TRUE"="Pass (>=80%)","FALSE"="Fail (<80%)"),
                    name="QC Status") +
  scale_y_continuous(limits=c(0,105), expand=expansion(mult=c(0,0.02)),
                     labels=function(x) paste0(x,"%")) +
  labs(title=expression(italic("K. pneumoniae")~"Kraken2 Taxonomic Classification QC"),
       subtitle=paste0("n = ", nrow(kraken), " isolates; dashed line = 80% classification threshold"),
       x="Isolate", y=expression(italic("K. pneumoniae")~"Reads (%)")) +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        legend.position=c(0.15, 0.25),
        legend.background=element_rect(colour="grey85", linewidth=0.3, fill="white"),
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"FigS1_Kraken2_QC.png"), figS1, width=12, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"FigS1_Kraken2_QC.pdf"), figS1, width=12, height=6)
message("FigS1 saved")
