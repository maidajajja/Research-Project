suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(forcats)
  library(scales)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pub_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, colour = "grey40", margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 12),
    axis.text        = element_text(size = 11, colour = "grey20"),
    axis.line        = element_line(colour = "grey30", linewidth = 0.5),
    axis.ticks       = element_line(colour = "grey30", linewidth = 0.4),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(12, 16, 10, 12)
  )

# ---- READ KLEBORATE ----
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
message("Kleborate columns: ", paste(colnames(kleb)[1:10], collapse=", "))

st_col  <- grep("mlst__ST",                        colnames(kleb), value=TRUE)[1]
vir_col <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]

kleb$ST  <- gsub("-.*", "", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))

# ---- FIG 2: ST DISTRIBUTION ----
important_hvkp <- c("ST23","ST65","ST380")
important_mdr  <- c("ST11","ST258","ST86","ST29","ST512","ST76")

st_counts <- kleb %>% count(ST, name="n") %>% arrange(desc(n))
top_sts   <- head(st_counts$ST, 20)
other_n   <- sum(st_counts$n[!st_counts$ST %in% top_sts])

st_df <- st_counts %>%
  filter(ST %in% top_sts) %>%
  bind_rows(data.frame(ST="Other STs", n=other_n)) %>%
  mutate(
    ST = factor(ST, levels=c(top_sts,"Other STs")),
    colour_group = case_when(
      ST %in% important_hvkp ~ "hvKp-associated",
      ST %in% important_mdr  ~ "MDR-associated",
      ST == "Other STs"      ~ "Other (grouped)",
      TRUE                   ~ "Other named ST"
    )
  )

colour_map <- c("hvKp-associated"="#C0392B","MDR-associated"="#2471A3",
                "Other named ST"="#5D6D7E","Other (grouped)"="#BDC3C7")

fig2 <- ggplot(st_df, aes(x=ST, y=n, fill=colour_group)) +
  geom_col(width=0.72, colour="white", linewidth=0.3) +
  geom_text(aes(label=n), vjust=-0.45, size=3.2, fontface="bold", colour="grey20") +
  scale_fill_manual(values=colour_map, name="ST category") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(title=expression("Sequence Type Distribution of 234 "*italic("K. pneumoniae")*" Isolates"),
       x="Sequence Type", y="Number of Isolates (n)") +
  pub_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9.5),
        legend.position=c(0.82,0.78),
        legend.background=element_rect(colour="grey85", linewidth=0.3, fill="white"))

ggsave(file.path(OUT,"Fig2_ST_distribution.png"), fig2, width=11, height=6, dpi=300, bg="white")
message("Fig2 saved")

# ---- FIG 3A: VIRULENCE DISTRIBUTION ----
vir_pal <- c("#F5EEF8","#C39BD3","#9B59B6","#7D3C98","#5B2C6F","#3B1A5A")
vir_df <- kleb %>% filter(!is.na(Vir)) %>% count(Vir, name="n") %>%
  mutate(pct=round(n/sum(n)*100,1), label=paste0(n,"\n(",pct,"%)"), Vir_f=factor(Vir))

fig3a <- ggplot(vir_df, aes(x=Vir_f, y=n, fill=Vir_f)) +
  geom_col(width=0.65, colour="white", linewidth=0.3) +
  geom_text(aes(label=label), vjust=-0.3, size=3.3, fontface="bold", colour="grey20", lineheight=0.95) +
  scale_fill_manual(values=vir_pal, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Virulence Score Distribution",
       subtitle=expression(italic("K. pneumoniae")*" (n = 234); Kleborate v3"),
       x="Virulence Score (0-5)", y="Number of Isolates (n)") +
  pub_theme

ggsave(file.path(OUT,"Fig3a_virulence_distribution.png"), fig3a, width=7, height=5.5, dpi=300, bg="white")
message("Fig3a saved")

# ---- FIG 3B: RESISTANCE DISTRIBUTION ----
res_pal <- c("#FDFEFE","#F1948A","#E74C3C","#C0392B")
res_df <- kleb %>% filter(!is.na(Res)) %>% count(Res, name="n") %>%
  mutate(pct=round(n/sum(n)*100,1), label=paste0(n,"\n(",pct,"%)"), Res_f=factor(Res))

fig3b <- ggplot(res_df, aes(x=Res_f, y=n, fill=Res_f)) +
  geom_col(width=0.55, colour="white", linewidth=0.3) +
  geom_text(aes(label=label), vjust=-0.3, size=3.3, fontface="bold", colour="grey20", lineheight=0.95) +
  scale_fill_manual(values=res_pal, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Antimicrobial Resistance Score Distribution",
       subtitle=expression(italic("K. pneumoniae")*" (n = 234); Kleborate v3"),
       x="Resistance Score (0-3)", y="Number of Isolates (n)") +
  pub_theme

ggsave(file.path(OUT,"Fig3b_resistance_distribution.png"), fig3b, width=7, height=5.5, dpi=300, bg="white")
message("Fig3b saved")

# ---- FIG 4: VIRULENCE VS RESISTANCE SCATTER ----
scatter_df <- kleb %>% filter(!is.na(Vir), !is.na(Res)) %>%
  mutate(Convergent=Vir>=4 & Res>=2,
         group=if_else(Convergent,"Convergent (n=20)","Non-convergent"))

bubble_df <- scatter_df %>% count(Res, Vir, group, name="count")

fig4 <- ggplot() +
  annotate("rect", xmin=1.7, xmax=3.3, ymin=3.7, ymax=5.3,
           fill="#FDEDEC", colour="#C0392B", linetype="dashed", linewidth=0.7, alpha=0.6) +
  annotate("text", x=2.5, y=5.25, label="Convergent zone (n = 20)",
           colour="#C0392B", fontface="bold", size=3.6, vjust=0) +
  geom_point(data=filter(bubble_df, group=="Non-convergent"),
             aes(x=Res, y=Vir, size=count),
             shape=21, fill="#AED6F1", colour="#2471A3", alpha=0.8, stroke=0.5) +
  geom_point(data=filter(bubble_df, group!="Non-convergent"),
             aes(x=Res, y=Vir, size=count),
             shape=24, fill="#E74C3C", colour="#922B21", alpha=0.9, stroke=0.6) +
  scale_size_area(max_size=12, name="No. isolates", breaks=c(1,5,10,20,40)) +
  scale_x_continuous(breaks=0:3, limits=c(-0.4,3.6)) +
  scale_y_continuous(breaks=0:5, limits=c(-0.5,5.6)) +
  labs(title=expression("Virulence vs. Resistance Scores in 234 "*italic("K. pneumoniae")*" Isolates"),
       subtitle="Bubble size proportional to number of isolates at each coordinate",
       x="Resistance Score (0-3)", y="Virulence Score (0-5)") +
  pub_theme + theme(legend.position="right")

ggsave(file.path(OUT,"Fig4_virulence_vs_resistance.png"), fig4, width=8.5, height=7, dpi=300, bg="white")
message("Fig4 saved")

# ---- FIG 5: PAN-GENOME ----
pg_df <- data.frame(
  Category=factor(c("Core\n(>=99% isolates)","Soft Core\n(95-99%)","Shell\n(15-95%)","Cloud\n(<15%)"),
                  levels=c("Core\n(>=99% isolates)","Soft Core\n(95-99%)","Shell\n(15-95%)","Cloud\n(<15%)")),
  Genes=c(3154,918,1604,18023),
  fill=c("#1A5276","#2E86C1","#E67E22","#C0392B")
)
pg_df$pct   <- round(pg_df$Genes/sum(pg_df$Genes)*100,1)
pg_df$label <- paste0(format(pg_df$Genes,big.mark=","),"\n(",pg_df$pct,"%)")

fig5 <- ggplot(pg_df, aes(x=Category, y=Genes, fill=fill)) +
  geom_col(width=0.6, colour="white", linewidth=0.3) +
  geom_text(aes(label=label), vjust=-0.3, size=3.4, fontface="bold", colour="grey20", lineheight=0.95) +
  scale_fill_identity() +
  scale_y_continuous(expand=expansion(mult=c(0,0.16)), labels=scales::comma) +
  labs(title=expression("Pan-genome Composition of 234 "*italic("K. pneumoniae")*" Isolates"),
       subtitle="Total pan-genome: 23,699 genes (Panaroo v1.3.4)",
       x="Gene Category", y="Number of Genes") +
  pub_theme

ggsave(file.path(OUT,"Fig5_pangenome.png"), fig5, width=8, height=6, dpi=300, bg="white")
message("Fig5 saved")

# ---- FIG 6: AMR PREVALENCE ----
amr <- read.table("/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder/amrfinder_all.tsv",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE)
message("AMR rows: ", nrow(amr))
message("AMR cols: ", paste(colnames(amr), collapse=" | "))

n_genomes <- 234
amr$Gene      <- amr[["Element symbol"]]
amr$DrugClass <- amr[["Class"]]

amr_prev <- amr %>%
  filter(!is.na(Gene), Gene != "") %>%
  count(Gene, DrugClass, name="n_hits") %>%
  mutate(prevalence=pmin(n_hits/n_genomes*100, 100)) %>%
  arrange(desc(prevalence)) %>%
  slice_head(n=20) %>%
  mutate(Gene=fct_reorder(Gene, prevalence))

class_colours <- c(
  "BETA-LACTAM"="#C0392B","AMINOGLYCOSIDE"="#2471A3",
  "TETRACYCLINE"="#D4AC0D","SULFONAMIDE"="#27AE60",
  "TRIMETHOPRIM"="#1A8A4A","QUINOLONE"="#8E44AD",
  "PHENICOL"="#E67E22","FOSFOMYCIN"="#F39C12","Other"="#7F8C8D"
)
amr_prev$DrugClass[!amr_prev$DrugClass %in% names(class_colours)] <- "Other"

fig6 <- ggplot(amr_prev, aes(x=prevalence, y=Gene, fill=DrugClass)) +
  geom_col(width=0.7, colour="white", linewidth=0.3) +
  geom_text(aes(label=paste0(n_hits," (",round(prevalence,1),"%)")),
            hjust=-0.08, size=3.1, colour="grey20") +
  scale_fill_manual(values=class_colours, name="Drug class") +
  scale_x_continuous(limits=c(0,80), expand=expansion(mult=c(0,0.25)),
                     labels=function(x) paste0(x,"%")) +
  labs(title=expression("Prevalence of AMR Genes in 234 "*italic("K. pneumoniae")*" Isolates"),
       subtitle="Top 20 genes by prevalence (AMRFinderPlus v4.2.7)",
       x="Prevalence (% of isolates)", y="AMR Gene") +
  pub_theme +
  theme(axis.text.y=element_text(face="italic", size=10),
        legend.position=c(0.78,0.28),
        legend.background=element_rect(colour="grey85", linewidth=0.3, fill="white"))

ggsave(file.path(OUT,"Fig6_AMR_prevalence.png"), fig6, width=9, height=8, dpi=300, bg="white")
message("Fig6 saved")

message("\n=== All figures complete ===")
