suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# NOTE: switched from the standalone mlst tool output to Kleborate's ST
# column. The standalone mlst tool was only ever run on 82/229 genomes
# (the BV-BRC/SRA-sourced ones) - it was never run on the 147 GCA/GCF
# genomes sourced from NCBI. Kleborate covers all genomes and its ST
# calls were already cross-validated against the standalone mlst tool
# (190 comparable genomes, 0 genuine disagreements - see Methods 2.4).
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
kleb$ST[kleb$ST==""] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

# Build a Sample ID that actually matches Kleborate's strain column.
# Same Assembly/SRA Accession fallback validated for Fig2/Fig3.
meta$GenomeID_str <- as.character(meta[["Genome ID"]])
meta$Sample <- ifelse(meta$GenomeID_str %in% kleb$strain, meta$GenomeID_str, NA)
unmatched_idx <- which(is.na(meta$Sample))
for (i in unmatched_idx) {
  acc <- meta[["Assembly Accession"]][i]
  sra <- meta[["SRA Accession"]][i]
  if (!is.na(acc) && acc != "") {
    m <- kleb$strain[startsWith(kleb$strain, acc)]
    if (length(m) >= 1) meta$Sample[i] <- m[1]
  } else if (!is.na(sra) && sra != "") {
    m <- kleb$strain[grepl(sra, kleb$strain, fixed = TRUE)]
    if (length(m) >= 1) meta$Sample[i] <- m[1]
  }
}

meta$Year    <- as.numeric(meta[["Collection Year"]])
meta$Country <- meta[["Isolation Country"]]
meta <- meta[!is.na(meta$Sample), c("Sample","Year","Country")]

st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, stringsAsFactors=FALSE)

df <- merge(st_map, meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$ST) & !is.na(df$Year) & df$Year > 2000, ]
df <- df[!is.na(df$Country) & df$Country != "" & df$Country != "NA" & df$Country != "Unknown", ]

top_sts <- names(sort(table(df$ST), decreasing=TRUE))[1:10]
st_labels <- paste0("ST", top_sts)

cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#000000","#888888","#44AA99")
st_cols <- setNames(cbf_pal[1:length(st_labels)], st_labels)

pub_theme <- theme_classic(base_size=14) +
  theme(axis.title=element_text(face="bold", size=13),
        axis.text=element_text(size=11, colour="grey20"),
        legend.title=element_text(face="bold", size=11),
        legend.text=element_text(size=10),
        plot.margin=margin(12,16,12,12))

# ── Panel A: Bubble plot ─────────────────────────────────────────────────────
st_year <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = factor(paste0("ST", ST), levels=rev(st_labels))) %>%
  group_by(Year, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  filter(n > 0)

figA <- ggplot(st_year, aes(x=Year, y=ST_label, size=n, fill=ST_label)) +
  geom_point(shape=21, colour="grey30", stroke=0.4, alpha=0.9) +
  geom_text(aes(label=n), colour="white", size=3, fontface="bold") +
  scale_size_continuous(range=c(6, 18), name="No. isolates") +
  scale_fill_manual(values=st_cols, guide="none") +
  scale_x_continuous(breaks=seq(min(st_year$Year), max(st_year$Year), by=2)) +
  labs(x="Collection Year", y="Sequence Type") +
  pub_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_line(colour="grey92", linewidth=0.4))

# ── Panel B: Stacked bar with total labels ───────────────────────────────────
df_cty <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = factor(paste0("ST", ST), levels=st_labels))

top_countries <- names(sort(table(df_cty$Country), decreasing=TRUE))[1:8]
df_cty <- df_cty %>% filter(Country %in% top_countries)
df_cty$Country <- factor(df_cty$Country, levels=top_countries)

totals <- df_cty %>% group_by(Country) %>% summarise(total=n(), .groups="drop")

figB <- ggplot(df_cty, aes(x=Country, fill=ST_label)) +
  geom_bar(position="stack", width=0.65, colour="white", linewidth=0.3) +
  geom_text(data=totals, aes(x=Country, y=total, label=total),
            inherit.aes=FALSE, vjust=-0.4, fontface="bold", size=3.5) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
  labs(x="Isolation Country", y="No. isolates") +
  pub_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4))

ggsave(file.path(OUT,"Fig5a_ST_over_time.png"),  figA, width=10, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5b_ST_by_country.png"), figB, width=8,  height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time.pdf"),  figA, width=10, height=6)
ggsave(file.path(OUT,"Fig5b_ST_by_country.pdf"), figB, width=8,  height=6)
message("Fig5 v6 saved - switched to Kleborate ST with Assembly/SRA fallback join")
