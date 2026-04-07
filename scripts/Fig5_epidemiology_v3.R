suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

mlst <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST/mlst_results.tsv",
                   sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(mlst) <- c("File","Scheme","ST","a1","a2","a3","a4","a5","a6","a7")
mlst$Sample <- gsub(".fasta$|.fna$|.fa$","", mlst$File)
mlst$ST <- as.character(mlst$ST)
mlst$ST[mlst$ST=="-"] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Year    <- as.numeric(meta[["Collection Year"]])
meta$Country <- meta[["Isolation Country"]]
meta <- meta[, c("Sample","Year","Country")]

df <- merge(mlst[,c("Sample","ST")], meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$ST) & !is.na(df$Year) & df$Year > 2000, ]
df <- df[!is.na(df$Country) & df$Country != "" & df$Country != "NA", ]

top_sts <- names(sort(table(df$ST), decreasing=TRUE))[1:10]
df$ST_group <- ifelse(df$ST %in% top_sts, paste0("ST",df$ST), "Other")
df$ST_group <- factor(df$ST_group, levels=c(paste0("ST",top_sts), "Other"))

st_pal <- c(brewer.pal(9,"Set1"), "#FF7F00", "grey80")
st_cols <- setNames(st_pal[1:length(levels(df$ST_group))], levels(df$ST_group))

pub_theme <- theme_classic(base_size=14) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        plot.subtitle=element_text(size=11, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=13),
        axis.text=element_text(size=11, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        legend.title=element_text(face="bold", size=11),
        legend.text=element_text(size=10),
        plot.margin=margin(12,16,12,12))

# Panel A: ST over time
figA <- ggplot(df, aes(x=Year, fill=ST_group)) +
  geom_bar(width=0.8, colour="white", linewidth=0.2) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_x_continuous(breaks=seq(min(df$Year),max(df$Year),by=2)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title=expression(italic("K. pneumoniae")~"Sequence Type Distribution Over Time"),
       subtitle="n = 234 isolates; collection year from BV-BRC metadata",
       x="Collection Year", y="Number of Isolates") +
  pub_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")

# Panel B: ST x Country heatmap - remove NA/Unknown
df_clean <- df[!is.na(df$Country) & df$Country != "NA" & df$Country != "" & df$Country != "Unknown", ]
top_countries <- names(sort(table(df_clean$Country), decreasing=TRUE))[1:8]
df_cty <- df_clean[df_clean$Country %in% top_countries & df_clean$ST %in% top_sts, ]
df_cty$ST_label <- paste0("ST", df_cty$ST)

st_country <- df_cty %>%
  group_by(Country, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  complete(Country=top_countries,
           ST_label=paste0("ST",top_sts),
           fill=list(n=0))
st_country$Country <- factor(st_country$Country, levels=rev(top_countries))
st_country$ST_label <- factor(st_country$ST_label, levels=paste0("ST",top_sts))

figB <- ggplot(st_country, aes(x=ST_label, y=Country, fill=n)) +
  geom_tile(colour="white", linewidth=0.8) +
  geom_text(aes(label=ifelse(n>0,n,"")), size=3.5, fontface="bold", colour="white") +
  scale_fill_gradientn(
    colours=c("#440154","#3B528B","#21908C","#5DC863","#FDE725"),
    name="No. isolates") +
  labs(title=expression(italic("K. pneumoniae")~"ST Distribution by Isolation Country"),
       subtitle="Top 10 STs across top 8 countries; n = 234 isolates",
       x="Sequence Type", y="Isolation Country") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(face="bold", size=14, hjust=0.5),
        plot.subtitle=element_text(size=11, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=13),
        axis.text=element_text(size=11, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(face="bold", size=11),
        plot.margin=margin(12,16,12,12))

ggsave(file.path(OUT,"Fig5a_ST_over_time.png"), figA, width=11, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5b_ST_by_country.png"), figB, width=11, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time.pdf"), figA, width=11, height=6)
ggsave(file.path(OUT,"Fig5b_ST_by_country.pdf"), figB, width=11, height=6)
message("Fig5 v2 saved")
