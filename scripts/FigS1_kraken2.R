suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"

df <- read.table("/scratch/users/k22017808/KP_Research_Project/plots/kraken2_summary.tsv",
                 sep="\t", header=TRUE, stringsAsFactors=FALSE)

df <- df %>% arrange(desc(KP_pct)) %>%
  mutate(Sample_idx = row_number(),
         Pass = KP_pct >= 50)

n_pass <- sum(df$Pass)
n_fail <- sum(!df$Pass)

fig <- ggplot(df, aes(x=Sample_idx, y=KP_pct, colour=Pass)) +
  geom_hline(yintercept=50, linetype="dashed", colour="#D55E00", linewidth=0.8) +
  geom_point(size=1.5, alpha=0.8) +
  scale_colour_manual(values=c("TRUE"="#0072B2","FALSE"="#D55E00"),
                      labels=c("TRUE"=paste0("Pass (n=",n_pass,")"),"FALSE"=paste0("Fail (n=",n_fail,")")),
                      name="QC status") +
  scale_y_continuous(limits=c(0,101), breaks=seq(0,100,20),
                     labels=function(x) paste0(x,"%")) +
  scale_x_continuous(expand=expansion(mult=c(0.01,0.01))) +
  labs(x="Sample (ordered by K. pneumoniae %)",
       y=expression(italic("K. pneumoniae")~"classified (%)")) +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(face="bold", size=10),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"FigS1_kraken2_QC.png"), fig,
       width=10, height=5, dpi=600, bg="white")
ggsave(file.path(OUT,"FigS1_kraken2_QC.pdf"), fig,
       width=10, height=5)
message("FigS1 saved at 600 DPI")
