suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load MLST results
mlst <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST/mlst_results.tsv",
                   sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(mlst) <- c("File","Scheme","ST","allele1","allele2","allele3","allele4","allele5","allele6","allele7")
mlst$Sample <- gsub(".fasta$|.fna$|.fa$","", mlst$File)
mlst$ST <- as.character(mlst$ST)
mlst$ST[mlst$ST=="-"] <- "Unknown"

# Load metadata
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Year    <- as.numeric(meta[["Collection Year"]])
meta$Country <- meta[["Isolation Country"]]
meta$Country[is.na(meta$Country)|meta$Country==""] <- "Unknown"
meta <- meta[, c("Sample","Year","Country")]

# Merge
df <- merge(mlst[,c("Sample","ST")], meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$Year) & df$Year > 2000 & df$ST != "Unknown", ]

# Top STs
top_sts <- names(sort(table(df$ST), decreasing=TRUE))[1:10]
df$ST_group <- ifelse(df$ST %in% top_sts, df$ST, "Other")
df$ST_group <- factor(df$ST_group, levels=c(top_sts, "Other"))

st_pal <- c(brewer.pal(9,"Set1"), "#FF7F00", "grey80")
st_cols <- setNames(st_pal[1:length(levels(df$ST_group))], levels(df$ST_group))

pub_theme <- theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

# Panel A: ST distribution over time
figA <- ggplot(df, aes(x=Year, fill=ST_group)) +
  geom_bar(width=0.8, colour="white", linewidth=0.2) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_x_continuous(breaks=seq(min(df$Year),max(df$Year),by=2)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title=expression(italic("K. pneumoniae")~"ST Distribution Over Time"),
       subtitle="Collection year from BV-BRC metadata",
       x="Collection Year", y="Number of Isolates") +
  pub_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")

# Panel B: ST x Country heatmap
top_countries <- names(sort(table(df$Country[df$Country!="Unknown"]), decreasing=TRUE))[1:8]
df_cty <- df[df$Country %in% top_countries & df$ST %in% top_sts, ]

st_country <- df_cty %>%
  group_by(Country, ST_group) %>%
  summarise(n=n(), .groups="drop") %>%
  complete(Country=top_countries, ST_group=factor(top_sts, levels=top_sts), fill=list(n=0))
st_country$Country <- factor(st_country$Country, levels=rev(top_countries))

figB <- ggplot(st_country, aes(x=ST_group, y=Country, fill=n)) +
  geom_tile(colour="white", linewidth=0.5) +
  scale_fill_gradientn(
    colours=c("#440154","#3B528B","#21908C","#5DC863","#FDE725"),
    name="No. isolates") +
  labs(title=expression(italic("K. pneumoniae")~"ST Distribution by Country"),
       subtitle="Top 10 STs x top 8 countries",
       x="Sequence Type", y="Isolation Country") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig5a_ST_over_time.png"), figA, width=11, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5b_ST_by_country.png"), figB, width=10, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time.pdf"), figA, width=11, height=6)
ggsave(file.path(OUT,"Fig5b_ST_by_country.pdf"), figB, width=10, height=6)
message("Fig5 saved")
