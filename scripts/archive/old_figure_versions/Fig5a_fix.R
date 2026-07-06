suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"

# NOTE: switched from the standalone mlst tool output to Kleborate's ST
# column - see Fig5_epidemiology_v7_klebfix.R for full rationale. The
# standalone mlst tool only covered 82/229 genomes.
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
kleb$ST <- sub("^ST", "", kleb$ST)
kleb$ST[kleb$ST==""] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

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
df <- df[!is.na(df$Country) & df$Country != "" &
         df$Country != "NA" & df$Country != "Unknown", ]

top_sts <- names(sort(table(df$ST), decreasing=TRUE))[1:10]
st_labels <- paste0("ST", top_sts)

cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#4D4D4D","#888888","#44AA99")
st_cols <- setNames(cbf_pal[1:length(st_labels)], st_labels)

st_year <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = factor(paste0("ST", ST), levels=rev(st_labels))) %>%
  group_by(Year, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  filter(n > 0)

# IMPORTANT: the original hardcoded position vector (c(1,2,3,5,6.5) for
# 2007/2009/2017/2021/2022) assumed exactly 5 distinct years from the
# old, incomplete 82-genome dataset. With full Kleborate coverage the
# set of distinct years present may differ, so positions are now built
# dynamically from whatever years are actually in the data, preserving
# the "extra gap before the most recent years" visual style as closely
# as possible rather than assuming a fixed 5-year set.
year_map <- sort(unique(st_year$Year))
n_years <- length(year_map)
if (n_years >= 2) {
  pos <- seq(1, n_years + 1.5, length.out = n_years)
  # widen the gap before the final year, mirroring the original intent
  pos[n_years] <- pos[n_years] + 0.5
} else {
  pos <- seq_len(n_years)
}
names(pos) <- as.character(year_map)
st_year$x_pos <- pos[as.character(st_year$Year)]

st_year$text_col <- ifelse(st_year$ST_label == "ST198", "grey20", "white")

figA <- ggplot(st_year, aes(x=x_pos, y=ST_label, size=n, fill=ST_label)) +
  geom_point(shape=21, colour="grey40", stroke=0.5, alpha=0.92) +
  geom_text(aes(label=n, colour=text_col), size=3.2, fontface="bold") +
  scale_size_continuous(range=c(7, 20), name="No. isolates",
                        breaks=c(1,3,6,9,12)) +
  scale_fill_manual(values=st_cols, guide="none") +
  scale_colour_identity() +
  scale_x_continuous(breaks=pos, labels=names(pos),
                     expand=expansion(add=c(0.4, 0.8))) +
  labs(x="Collection Year", y="Sequence Type") +
  theme_classic(base_size=14) +
  theme(axis.title=element_text(face="bold", size=13),
        axis.text=element_text(size=11, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_line(colour="grey85", linewidth=0.5),
        panel.grid.minor=element_blank(),
        legend.title=element_text(face="bold", size=11),
        legend.text=element_text(size=10),
        legend.background=element_blank(),
        legend.key=element_rect(fill="grey92", colour=NA),
        plot.margin=margin(12,24,12,12))

ggsave(file.path(OUT,"Fig5a_ST_over_time.png"), figA,
       width=9, height=6.5, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time.pdf"), figA,
       width=9, height=6.5)
message("Fig5a fix3 saved - switched to Kleborate ST, dynamic year positions")
