suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(forcats)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv", stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample <- as.character(meta[["Genome ID"]])
meta$Year <- as.numeric(meta[["Collection Year"]])
meta$Country <- meta[["Isolation Country"]]
meta$Country[is.na(meta$Country)|meta$Country==""] <- "Unknown"
meta <- meta[!is.na(meta$Year) & meta$Year > 2000, ]

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, stringsAsFactors=FALSE)

df <- merge(meta, st_map, by="Sample", all.x=TRUE)
df$ST[is.na(df$ST)] <- "Unknown"

top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE))[1:7]
df$ST_group <- ifelse(df$ST %in% top_sts, df$ST, "Other")
df$ST_group <- factor(df$ST_group, levels=c(top_sts, "Other"))

st_pal <- brewer.pal(7,"Set1")
st_cols <- setNames(c(st_pal,"grey80"), c(top_sts,"Other"))

# Panel A: ST distribution over time
fig5a <- ggplot(df, aes(x=Year, fill=ST_group)) +
  geom_bar(width=0.8, colour="white", linewidth=0.2) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_x_continuous(breaks=seq(min(df$Year, na.rm=TRUE), max(df$Year, na.rm=TRUE), by=2)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title=expression(italic("K. pneumoniae")~"Sequence Type Distribution Over Time"),
       x="Collection Year", y="Number of Isolates") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        legend.position="right",
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

# Panel B: ST distribution by country
country_counts <- df %>% count(Country) %>% arrange(desc(n))
top_countries <- head(country_counts$Country, 8)
df$Country_group <- ifelse(df$Country %in% top_countries, df$Country, "Other")
df$Country_group <- factor(df$Country_group, levels=c(top_countries, "Other"))

fig5b <- ggplot(df, aes(x=fct_infreq(Country_group), fill=ST_group)) +
  geom_bar(width=0.7, colour="white", linewidth=0.2) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title=expression(italic("K. pneumoniae")~"Sequence Type Distribution by Country"),
       x="Isolation Country", y="Number of Isolates") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major.y=element_line(colour="grey92", linewidth=0.4),
        legend.position="right",
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig5a_ST_by_year.png"), fig5a, width=10, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5b_ST_by_country.png"), fig5b, width=10, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_by_year.pdf"), fig5a, width=10, height=6)
ggsave(file.path(OUT,"Fig5b_ST_by_country.pdf"), fig5b, width=10, height=6)
message("Fig5 saved")
