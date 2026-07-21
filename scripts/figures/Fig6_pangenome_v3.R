suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)

})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

rtab_path <- "/scratch/users/k22017808/KP_Research_Project/06_Pangenome_final229/gene_presence_absence.Rtab"
rtab <- read.table(rtab_path, sep="\t", header=TRUE,
                   stringsAsFactors=FALSE, check.names=FALSE)
rownames(rtab) <- rtab[,1]
rtab <- rtab[,-1]
mat <- as.matrix(rtab)

n_samples <- ncol(mat)
n_genes   <- nrow(mat)
gene_freq <- rowSums(mat) / n_samples * 100

cat("Total genes:", n_genes, "\n")
cat("Total isolates:", n_samples, "\n")

# Pan-genome categories
core     <- sum(gene_freq >= 99)
softcore <- sum(gene_freq >= 95 & gene_freq < 99)
shell    <- sum(gene_freq >= 15 & gene_freq < 95)
cloud    <- sum(gene_freq < 15)
cat("Core:", core, "Soft core:", softcore, "Shell:", shell, "Cloud:", cloud, "\n")

# ── Fig6a: Pan-genome composition bar chart ──────────────────────────────────
pg <- data.frame(
  Category = factor(c("Core\n(>=99%)", "Soft core\n(95-99%)",
                       "Shell\n(15-95%)", "Cloud\n(<15%)"),
                    levels=c("Core\n(>=99%)","Soft core\n(95-99%)",
                             "Shell\n(15-95%)","Cloud\n(<15%)")),
  Genes = c(core, softcore, shell, cloud),
  Fill  = c("#1A5276","#2E86C1","#E67E22","#C0392B"))
pg$Pct   <- round(pg$Genes / sum(pg$Genes) * 100, 1)
pg$Label <- paste0(formatC(pg$Genes, format="d", big.mark=","),
                   "\n(", pg$Pct, "%)")

fig6a <- ggplot(pg, aes(x=Category, y=Genes, fill=Fill)) +
  geom_col(width=0.6, show.legend=FALSE) +
  geom_text(aes(label=Label), vjust=-0.3, size=4,
            fontface="bold", colour="grey20", lineheight=0.9) +
  scale_fill_identity() +
  scale_y_continuous(expand=expansion(mult=c(0,0.18)),
                     labels=scales::comma) +
  labs(x="Pan-genome category", y="Number of genes",
       subtitle=paste0("Total pan-genome: ",
                       formatC(n_genes, format="d", big.mark=","),
                       " genes | n=", n_samples, " isolates")) +
  theme_classic(base_size=12) +
  theme(
    axis.title = element_text(face="bold", size=11),
    axis.text  = element_text(size=10, colour="grey20"),
    plot.subtitle = element_text(size=9, colour="grey40"),
    panel.grid.major.y = element_line(colour="grey92", linewidth=0.4),
    plot.margin = margin(10,10,10,10))

# ── Fig6b: Gene frequency histogram with all 4 boundaries ────────────────────
freq_df <- data.frame(freq=gene_freq)

# Colour regions by category
fig6b <- ggplot(freq_df, aes(x=freq)) +
  # Coloured background regions
  annotate("rect", xmin=0,  xmax=15, ymin=0, ymax=Inf,
           fill="#C0392B", alpha=0.06) +
  annotate("rect", xmin=15, xmax=95, ymin=0, ymax=Inf,
           fill="#E67E22", alpha=0.06) +
  annotate("rect", xmin=95, xmax=99, ymin=0, ymax=Inf,
           fill="#2E86C1", alpha=0.08) +
  annotate("rect", xmin=99, xmax=100, ymin=0, ymax=Inf,
           fill="#1A5276", alpha=0.08) +
  geom_histogram(binwidth=2, fill="#2980B9", colour="white",
                 linewidth=0.2) +
  # Boundary lines - all 4 categories
  geom_vline(xintercept=15, linetype="dashed",
             colour="#C0392B", linewidth=0.8) +
  geom_vline(xintercept=95, linetype="dashed",
             colour="#E67E22", linewidth=0.8) +
  geom_vline(xintercept=99, linetype="dashed",
             colour="#1A5276", linewidth=0.8) +
  # Labels
  annotate("text", x=7,  y=Inf, label="Cloud\n(<15%)",
           colour="#C0392B", vjust=1.3, size=3.5, fontface="bold") +
  annotate("text", x=55, y=Inf, label="Shell\n(15-95%)",
           colour="#E67E22", vjust=1.3, size=3.5, fontface="bold") +
  annotate("text", x=92, y=3500, label="Soft
core",
           colour="#2E86C1", vjust=1.3, size=3, fontface="bold") +
  annotate("text", x=100, y=Inf, label="Core (>=99%)", colour="#1A5276", vjust=1.3, size=3, fontface="bold", hjust=1) +
  scale_x_continuous(breaks=seq(0,100,10),
                     labels=function(x) paste0(x,"%"),
                     limits=c(0,101)) +
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.05))) +
  labs(x="Gene frequency (% isolates)", y="Number of genes") +
  theme_classic(base_size=12) +
  theme(
    axis.title = element_text(face="bold", size=11),
    axis.text  = element_text(size=10, colour="grey20"),
    panel.grid.major.y = element_line(colour="grey92", linewidth=0.4),
    plot.margin = margin(10,10,10,10))

# ── Fig6c: Pan-genome accumulation curve ─────────────────────────────────────
# Rarefaction — random genome addition order, 50 permutations
set.seed(42)
n_perm <- 50
accum <- matrix(NA, nrow=n_perm, ncol=n_samples)
core_accum <- matrix(NA, nrow=n_perm, ncol=n_samples)

cat("Running", n_perm, "permutations for accumulation curve...\n")
for(p in 1:n_perm) {
  idx <- sample(n_samples)
  pan  <- 0
  core_genes <- 1:n_genes
  pan_vec  <- numeric(n_samples)
  core_vec <- numeric(n_samples)
  for(i in seq_along(idx)) {
    present <- which(mat[, idx[i]] == 1)
    pan <- union(pan, present)
    if(i == 1) {
      core_genes <- present
    } else {
      core_genes <- intersect(core_genes, present)
    }
    pan_vec[i]  <- length(pan)
    core_vec[i] <- length(core_genes)
  }
  accum[p,]      <- pan_vec
  core_accum[p,] <- core_vec
}

accum_df <- data.frame(
  n_genomes  = 1:n_samples,
  pan_mean   = colMeans(accum),
  pan_lo     = apply(accum, 2, quantile, 0.1),
  pan_hi     = apply(accum, 2, quantile, 0.9),
  core_mean  = colMeans(core_accum),
  core_lo    = apply(core_accum, 2, quantile, 0.1),
  core_hi    = apply(core_accum, 2, quantile, 0.9))

fig6c <- ggplot(accum_df, aes(x=n_genomes)) +
  # Pan-genome ribbon and line
  geom_ribbon(aes(ymin=pan_lo, ymax=pan_hi),
              fill="#C0392B", alpha=0.15) +
  geom_line(aes(y=pan_mean, colour="Pan-genome"),
            linewidth=1.2) +
  # Core genome ribbon and line
  geom_ribbon(aes(ymin=core_lo, ymax=core_hi),
              fill="#1A5276", alpha=0.15) +
  geom_line(aes(y=core_mean, colour="Core genome"),
            linewidth=1.2) +
  scale_colour_manual(
    values=c("Pan-genome"="#C0392B", "Core genome"="#1A5276"),
    name=NULL) +
  scale_y_continuous(labels=scales::comma,
                     expand=expansion(mult=c(0,0.05))) +
  scale_x_continuous(breaks=seq(0, n_samples, 25)) +
  labs(x="Number of genomes", y="Number of genes",
       subtitle=paste0("Pan-genome accumulation curve (", n_perm,
                       " permutations, shading = 10th-90th percentile)")) +
  theme_classic(base_size=12) +
  theme(
    axis.title = element_text(face="bold", size=11),
    axis.text  = element_text(size=10, colour="grey20"),
    plot.subtitle = element_text(size=9, colour="grey40"),
    legend.position = c(0.7, 0.5),
    legend.text = element_text(size=10),
    legend.background = element_blank(),
    panel.grid.major = element_line(colour="grey92", linewidth=0.4),
    plot.margin = margin(10,10,10,10))

# Save individual figures
ggsave(file.path(OUT,"Fig6a_pangenome_composition_v3.png"),
       fig6a, width=8, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6a_pangenome_composition_v3.pdf"),
       fig6a, width=8, height=6)

ggsave(file.path(OUT,"Fig6b_pangenome_frequency_v3.png"),
       fig6b, width=9, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6b_pangenome_frequency_v3.pdf"),
       fig6b, width=9, height=6)

ggsave(file.path(OUT,"Fig6c_pangenome_accumulation_v3.png"),
       fig6c, width=9, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6c_pangenome_accumulation_v3.pdf"),
       fig6c, width=9, height=6)

message("Fig6 v2 saved - composition, frequency, accumulation curve")

# ── Fig6ab: Combined panel ────────────────────────────────────────────────────
suppressPackageStartupMessages(library(patchwork))

fig6ab <- fig6a + fig6b +
  plot_layout(ncol=2, widths=c(1,1.4)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size=14, face="bold"))
  )

ggsave(file.path(OUT,"Fig6ab_pangenome_composition_frequency_v3.png"),
       fig6ab, width=16, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig6ab_pangenome_composition_frequency_v3.pdf"),
       fig6ab, width=16, height=7)
message("Fig6ab combined panel saved")
