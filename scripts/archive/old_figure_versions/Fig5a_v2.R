suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
kleb$ST <- sub("^ST", "", kleb$ST)
kleb$ST[kleb$ST==""] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

# Match samples
meta$GenomeID_str <- as.character(meta[["Genome ID"]])
meta$Sample <- ifelse(meta$GenomeID_str %in% kleb$strain, meta$GenomeID_str, NA)
unmatched_idx <- which(is.na(meta$Sample))
for(i in unmatched_idx) {
  acc <- meta[["Assembly Accession"]][i]
  sra <- meta[["SRA Accession"]][i]
  if(!is.na(acc) && acc != "") {
    m <- kleb$strain[startsWith(kleb$strain, acc)]
    if(length(m) >= 1) meta$Sample[i] <- m[1]
  } else if(!is.na(sra) && sra != "") {
    m <- kleb$strain[grepl(sra, kleb$strain, fixed=TRUE)]
    if(length(m) >= 1) meta$Sample[i] <- m[1]
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

# Top 7 STs only - consistent with other figures
top_sts <- c("512","29","65","258","86","11","23")
df$ST_label <- ifelse(df$ST %in% top_sts, paste0("ST", df$ST), "Other")
st_labels <- c(paste0("ST", top_sts), "Other")
df$ST_label <- factor(df$ST_label, levels=st_labels)

# Filter to top 7 only - remove Other for cleaner figure
df <- df %>% filter(ST_label != "Other")
df$ST_label <- droplevels(df$ST_label)
st_labels <- levels(df$ST_label)

st_year <- df %>%
  group_by(ST_label, Year) %>%
  summarise(n=n(), .groups="drop")

cat("Max isolates in one year-ST:", max(st_year$n), "\n")
cat("ST-year combinations:", nrow(st_year), "\n")

# Year axis positions
all_years <- seq(min(st_year$Year), max(st_year$Year))
pos <- setNames(all_years, all_years)

# Paul Tol muted palette - consistent with all other figures
st_cols <- setNames(c(
  "#77AADD",  # ST512 - medium blue
  "#AAAA00",  # ST29 - olive
  "#EE8866",  # ST65 - muted orange
  "#332288",  # ST258 - dark indigo
  "#CC99BB",  # ST86 - muted mauve
  "#225522",  # ST11 - dark forest green
  "#774411"   # ST23 - dark brown
), st_labels)

# Text colour - white on dark, dark on light
st_year$text_col <- ifelse(
  st_year$ST_label %in% c("ST23","ST11","ST258","ST29"), "white", "grey20")

max_n <- max(st_year$n)

fig5a <- ggplot(st_year, aes(x=Year, y=ST_label, size=n, fill=ST_label)) +
  geom_point(shape=21, colour="grey40", stroke=0.5, alpha=0.92) +
  geom_text(aes(label=n, colour=text_col), size=3.2, fontface="bold") +
  scale_size_area(max_size=18, name="No. isolates",
                  breaks=c(1, 5, 10, 20, max_n)) +
  scale_fill_manual(values=st_cols, guide="none") +
  scale_colour_identity() +
  scale_x_continuous(breaks=all_years,
                     labels=ifelse(all_years %% 2 == 1, all_years, ""),
                     limits=c(min(all_years)-0.5, max(all_years)+0.5)) +
  labs(x="Collection Year", y="Sequence Type") +
  theme_classic(base_size=13) +
  theme(
    axis.title = element_text(face="bold", size=12),
    axis.text  = element_text(size=10, colour="grey20"),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_line(colour="grey88", linewidth=0.4),
    legend.position = "right",
    legend.title = element_text(face="bold", size=10),
    legend.text  = element_text(size=9),
    legend.key   = element_rect(fill="grey92", colour=NA),
    plot.margin  = margin(10,15,10,10))

ggsave(file.path(OUT,"Fig5a_ST_over_time_v2.png"),
       fig5a, width=12, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time_v2.pdf"),
       fig5a, width=12, height=6)

message("Fig5a v2 saved")
