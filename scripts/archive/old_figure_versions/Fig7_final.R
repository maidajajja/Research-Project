suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample <- as.character(meta[["Genome ID"]])
meta$Source <- meta[["Isolation Source"]]
meta$Source[is.na(meta$Source)|meta$Source==""] <- "Unknown"
meta <- meta[, c("Sample","Source")]

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col  <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST  <- gsub("-.*","", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))

df <- merge(data.frame(Sample=kleb$strain, ST=kleb$ST,
                       Vir=kleb$Vir, Res=kleb$Res),
            meta, by="Sample", all.x=TRUE)
df$ST[is.na(df$ST)] <- "Unknown"
df <- df[!is.na(df$Vir) & !is.na(df$Res), ]
df$Convergent <- df$Vir >= 4 & df$Res >= 2

top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE))[1:7]
df$ST_group <- ifelse(df$ST %in% top_sts, paste0("ST", df$ST), "Other")
st_labels <- c(paste0("ST", top_sts), "Other")
df$ST_group <- factor(df$ST_group, levels=st_labels)

# CBF palette
cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","grey60")
st_cols <- setNames(cbf_pal[1:length(st_labels)], st_labels)

bubble_df <- df %>%
  group_by(Res, Vir, Convergent, ST_group) %>%
  summarise(count=n(), .groups="drop")

n_conv <- sum(df$Convergent)

fig7 <- ggplot() +
  # convergent zone shading
  annotate("rect", xmin=1.6, xmax=3.4, ymin=3.6, ymax=5.4,
           fill="#FEE8E7", colour="#CC0000", linetype="dashed",
           linewidth=0.7, alpha=0.6) +
  annotate("text", x=2.5, y=5.45,
           label=paste0("Convergent zone (n = ", n_conv, ")"),
           colour="#CC0000", fontface="bold", size=3.8, vjust=0) +
  # non-convergent circles
  geom_point(data=filter(bubble_df, !Convergent),
             aes(x=Res, y=Vir, size=count, fill=ST_group),
             shape=21, colour="grey40", alpha=0.85, stroke=0.5) +
  # convergent triangles
  geom_point(data=filter(bubble_df, Convergent),
             aes(x=Res, y=Vir, size=count, fill=ST_group),
             shape=24, colour="grey20", alpha=0.95, stroke=0.6) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_size_area(max_size=16, name="No. isolates",
                  breaks=c(1,5,10,20,40)) +
  scale_x_continuous(breaks=0:3, limits=c(-0.5,3.7),
                     labels=c("0","1","2","3")) +
  scale_y_continuous(breaks=0:5, limits=c(-0.7,5.8)) +
  # shape key as text annotation
  annotate("text", x=-0.45, y=-0.62,
           label="Circle = Non-convergent     Triangle = Convergent",
           size=3.3, colour="grey35", hjust=0, fontface="italic") +
  labs(x="Resistance score (0\u20133)", y="Virulence score (0\u20135)") +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major=element_line(colour="grey92", linewidth=0.4),
        legend.position="right",
        legend.title=element_text(face="bold", size=10),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        plot.margin=margin(12,15,12,12))

ggsave(file.path(OUT,"Fig7_convergence.png"), fig7,
       width=10, height=8, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig7_convergence.pdf"), fig7,
       width=10, height=8)
message("Fig7 final saved")
