suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(tidyr)
  library(RColorBrewer)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load plasmid data
plasmid <- read.table("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/plasmidfinder_results.tsv",
                      sep="\t", header=TRUE, stringsAsFactors=FALSE)
plasmid <- plasmid[plasmid$Plasmid != "None", ]

# Load metadata
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Country <- meta[["Isolation Country"]]
meta$Country[is.na(meta$Country)|meta$Country==""] <- "Unknown"
meta <- meta[, c("Sample","Country")]

# Load ST from Kleborate
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, stringsAsFactors=FALSE)

# Merge
df <- plasmid %>%
  left_join(st_map, by="Sample") %>%
  left_join(meta, by="Sample")
df$ST[is.na(df$ST)] <- "Unknown"
df$Country[is.na(df$Country)] <- "Unknown"

# Top STs
top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE))[1:7]
df$ST_group <- ifelse(df$ST %in% top_sts, df$ST, "Other")

# Top 15 plasmids
top_plasmids <- names(sort(table(df$Plasmid), decreasing=TRUE))[1:15]
df_filt <- df[df$Plasmid %in% top_plasmids, ]

# Count per plasmid per ST
plasmid_st <- df_filt %>%
  group_by(ST_group, Plasmid) %>%
  summarise(n=n_distinct(Sample), .groups="drop")
plasmid_st$Plasmid <- factor(plasmid_st$Plasmid,
  levels=names(sort(tapply(plasmid_st$n, plasmid_st$Plasmid, sum))))

# Count per plasmid per Country
top_countries <- names(sort(table(df_filt$Country[df_filt$Country != "Unknown" & !is.na(df_filt$Country)]), decreasing=TRUE))[1:6]
df_country <- df_filt[df_filt$Country %in% top_countries, ]
plasmid_country <- df_country %>%
  group_by(Country, Plasmid) %>%
  summarise(n=n_distinct(Sample), .groups="drop") %>%
  complete(Country, Plasmid=top_plasmids, fill=list(n=0))
plasmid_country$Plasmid <- factor(plasmid_country$Plasmid,
  levels=names(sort(tapply(plasmid_country$n, plasmid_country$Plasmid, sum))))

# ST colours
st_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
st_cols <- setNames(c(st_pal,"grey80"), c(top_sts,"Other"))

# Country colours
cty_cols <- setNames(colorRampPalette(brewer.pal(6,"Set2"))(length(top_countries)), top_countries)

pub_theme <- theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.y=element_text(face="italic"),
        panel.grid.major.x=element_line(colour="grey92", linewidth=0.4),
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

# Panel A: by ST
figA <- ggplot(plasmid_st, aes(x=n, y=Plasmid, fill=ST_group)) +
  geom_col(width=0.7, colour="white", linewidth=0.3) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) +
  labs(x="Number of Isolates", y="Plasmid Replicon Type") +
  pub_theme

# Panel B: by Country (heatmap style like supervisor)
figB <- ggplot(plasmid_country, aes(x=Plasmid, y=Country, fill=n)) +
  geom_tile(colour="white", linewidth=0.5) +
  scale_fill_gradientn(
    colours=c("#440154","#3B528B","#21908C","#5DC863","#FDE725"),
    name="No. isolates") +
  labs(x="Plasmid Replicon Type", y="Isolation Country") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.x=element_text(angle=45, hjust=1, face="italic"),
        legend.title=element_text(face="bold", size=10),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig4a_plasmid_by_ST.png"), figA, width=10, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4b_plasmid_by_country.png"), figB, width=10, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4a_plasmid_by_ST.pdf"), figA, width=10, height=7)
ggsave(file.path(OUT,"Fig4b_plasmid_by_country.pdf"), figB, width=10, height=6)
message("Fig4 saved")
