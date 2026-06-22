suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load Kleborate
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col  <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST  <- gsub("-.*","", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))

# Load metadata
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Source  <- meta[["Isolation Source"]]
meta$Health  <- meta[["Host Health"]]
meta$Source[is.na(meta$Source)|meta$Source==""] <- "Unknown"
meta$Health[is.na(meta$Health)|meta$Health==""] <- "Unknown"
meta <- meta[, c("Sample","Source","Health")]

df <- merge(data.frame(Sample=kleb$strain, ST=kleb$ST,
                       Vir=kleb$Vir, Res=kleb$Res),
            meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$Vir) & !is.na(df$Res), ]
df$Source[is.na(df$Source)] <- "Unknown"
df$Health[is.na(df$Health)] <- "Unknown"
df$ST[is.na(df$ST)] <- "Unknown"

# Define convergent
df$Convergent <- df$Vir >= 4 & df$Res >= 2

# Top STs
top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE))[1:7]
df$ST_group <- ifelse(df$ST %in% top_sts, paste0("ST",df$ST), "Other")

# Bubble counts
bubble_df <- df %>%
  group_by(Res, Vir, Convergent, ST_group, Source) %>%
  summarise(count=n(), .groups="drop")

st_pal <- brewer.pal(7,"Set1")
st_label_levels <- c(paste0("ST",top_sts),"Other")
st_cols <- setNames(c(st_pal[1:length(top_sts)],"grey70"), st_label_levels)

# Source shapes
sources <- unique(df$Source[df$Source != "Unknown"])
source_shapes <- setNames(c(21,22,23,24,25,21)[1:length(sources)], sources)

# Panel A: main convergence bubble plot coloured by ST
figA <- ggplot() +
  annotate("rect", xmin=1.7, xmax=3.3, ymin=3.7, ymax=5.3,
           fill="#FDEDEC", colour="#C0392B", linetype="dashed",
           linewidth=0.8, alpha=0.5) +
  annotate("text", x=2.5, y=5.25,
           label="Convergent zone (n = 20)",
           colour="#C0392B", fontface="bold", size=4, vjust=0) +
  geom_point(data=filter(bubble_df, !Convergent),
             aes(x=Res, y=Vir, size=count, fill=ST_group),
             shape=21, colour="grey40", alpha=0.85, stroke=0.5) +
  geom_point(data=filter(bubble_df, Convergent),
             aes(x=Res, y=Vir, size=count, fill=ST_group),
             shape=24, colour="#922B21", alpha=0.95, stroke=0.8) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_size_area(max_size=14, name="No. isolates", breaks=c(1,5,10,20,40)) +
  scale_x_continuous(breaks=0:3, limits=c(-0.5,3.6),
                     labels=c("0","1","2","3")) +
  scale_y_continuous(breaks=0:5, limits=c(-0.6,5.6)) +
  labs(title=expression(italic("K. pneumoniae")~"Virulence-Resistance Convergence (n = 234)"),
       subtitle="Convergent = virulence score >= 4 AND resistance score >= 2; triangle = convergent isolate",
       x="Resistance Score (0-3)", y="Virulence Score (0-5)") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major=element_line(colour="grey92", linewidth=0.4),
        legend.position="right",
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

# Panel B: convergent isolates breakdown by ST and Source
conv_df <- df[df$Convergent, ]
conv_df$ST_label <- paste0("ST", conv_df$ST)

top_sources <- names(sort(table(conv_df$Source), decreasing=TRUE))[1:5]
conv_df$Src_group <- ifelse(conv_df$Source %in% top_sources, conv_df$Source, "Other")

figB <- ggplot(conv_df, aes(x=ST_label, fill=Health)) +
  geom_bar(width=0.7, colour="white", linewidth=0.3) +
  scale_fill_brewer(palette="Set2", name="Host Health") +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), breaks=1:10) +
  labs(title=expression("Convergent "*italic("K. pneumoniae")*" Isolates by ST and Source (n = 20)"),
       subtitle="Coloured by host health status; virulence >= 4 and resistance >= 2",
       x="Sequence Type", y="Number of Isolates") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=11, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        legend.position="right",
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig7a_convergence_bubble.png"), figA, width=10, height=8, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig7b_convergent_breakdown.png"), figB, width=9, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig7a_convergence_bubble.pdf"), figA, width=10, height=8)
ggsave(file.path(OUT,"Fig7b_convergent_breakdown.pdf"), figB, width=9, height=6)
message("Fig7 saved")
