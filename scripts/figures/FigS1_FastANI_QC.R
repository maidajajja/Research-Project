suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

ani <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI/fastani_results.txt",
                  sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(ani) <- c("Query","Reference","ANI","Fragments","Total")

# Get average ANI per isolate (excluding self-comparisons)
clean_name <- function(x) {
  x <- basename(x)
  x <- gsub(".fna$","", x)
  x <- gsub("_genomic$","", x)
  x
}
ani$Query     <- clean_name(ani$Query)
ani$Reference <- clean_name(ani$Reference)

# Get mean ANI per isolate across all comparisons
ani_per_isolate <- ani %>%
  filter(Query != Reference) %>%
  group_by(Query) %>%
  summarise(mean_ANI=mean(ANI), min_ANI=min(ANI), max_ANI=max(ANI), .groups="drop")

message("Isolates: ", nrow(ani_per_isolate))
message("ANI range: ", round(min(ani_per_isolate$mean_ANI),3), " - ", round(max(ani_per_isolate$mean_ANI),3))

figS1 <- ggplot(ani_per_isolate, aes(x=mean_ANI)) +
  geom_histogram(binwidth=0.1, fill="#2471A3", colour="white", linewidth=0.2) +
  geom_vline(xintercept=95, linetype="dashed", colour="#C0392B", linewidth=0.8) +
  annotate("text", x=95.05, y=Inf, label="Species boundary (95% ANI)",
           colour="#C0392B", fontface="bold", size=3.5, hjust=0, vjust=2) +
  scale_x_continuous(labels=function(x) paste0(x,"%")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title=expression("FastANI Genomic Identity of 234 " * italic("K. pneumoniae") * " Isolates"),
       subtitle="Mean pairwise ANI against all other isolates; all isolates confirmed as K. pneumoniae (ANI > 95%)",
       x="Mean Average Nucleotide Identity (%)",
       y="Number of Isolates") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"FigS1_FastANI_QC.png"), figS1, width=10, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"FigS1_FastANI_QC.pdf"), figS1, width=10, height=6)
message("FigS1 saved")
