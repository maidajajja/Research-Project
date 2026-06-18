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
meta$Country <- meta[["Isolation Country"]]
meta <- meta[, c("Sample","Country")]
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, stringsAsFactors=FALSE)
plasmid <- read.table("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/plasmidfinder_results.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
plasmid <- plasmid[plasmid$Plasmid != "None", ]
df <- merge(plasmid, st_map, by="Sample", all.x=TRUE)
df <- merge(df, meta, by="Sample", all.x=TRUE)
df$ST[is.na(df$ST)] <- "Unknown"
df$Country[is.na(df$Country)] <- "Unknown"
top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE))[1:7]
df$ST_group <- ifelse(df$ST %in% top_sts, df$ST, "Other")
top_plasmids <- names(sort(table(df$Plasmid), decreasing=TRUE))[1:15]
df_filt <- df[df$Plasmid %in% top_plasmids, ]
plasmid_st <- aggregate(Sample ~ ST_group + Plasmid, data=df_filt, FUN=function(x) length(unique(x)))
colnames(plasmid_st)[3] <- "n_isolates"
plasmid_st$Plasmid <- factor(plasmid_st$Plasmid, levels=names(sort(tapply(plasmid_st$n_isolates, plasmid_st$Plasmid, sum))))
st_pal <- RColorBrewer::brewer.pal(7,"Set1")
st_cols <- setNames(c(st_pal,"grey80"), c(top_sts,"Other"))
fig4 <- ggplot(plasmid_st, aes(x=n_isolates, y=Plasmid, fill=ST_group)) +
  geom_col(width=0.7, colour="white", linewidth=0.3) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_x_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title=expression(italic("K. pneumoniae")~"Plasmid Replicon Distribution by Lineage"), subtitle="Top 15 plasmid types (PlasmidFinder v2.1.6)", x="Number of Isolates", y="Plasmid Replicon Type") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5), plot.subtitle=element_text(size=10, hjust=0.5, colour="grey40"), axis.title=element_text(face="bold", size=12), axis.text=element_text(size=10, colour="grey20"), axis.text.y=element_text(face="italic"), panel.grid.major.x=element_line(colour="grey92", linewidth=0.4), legend.position="right", legend.title=element_text(face="bold", size=10), plot.margin=margin(10,15,10,10))
ggsave(file.path(OUT,"Fig4_plasmid_distribution.png"), fig4, width=10, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4_plasmid_distribution.pdf"), fig4, width=10, height=7)
message("Fig4 saved")
