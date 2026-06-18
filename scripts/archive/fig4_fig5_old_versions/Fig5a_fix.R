suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"

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

actual_years <- sort(unique(st_year$Year))
st_year$text_col <- ifelse(st_year$ST_label == "ST198", "grey20", "white")

figA <- ggplot(st_year, aes(x=Year, y=ST_label, size=n, fill=ST_label)) +
  geom_point(shape=21, colour="grey40", stroke=0.5, alpha=0.92) +
  geom_text(aes(label=n, colour=text_col), size=3.2, fontface="bold") +
  scale_size_continuous(range=c(7, 20), name="No. isolates",
                        breaks=c(1,3,6,9,12)) +
  scale_fill_manual(values=st_cols, guide="none") +
  scale_colour_identity() +
  scale_x_continuous(breaks=actual_years, labels=actual_years,
                     expand=expansion(add=c(0.5, 1.2))) +
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
        legend.key=element_blank(),
        plot.margin=margin(12,24,12,12))

ggsave(file.path(OUT,"Fig5a_ST_over_time.png"), figA,
       width=9, height=6.5, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time.pdf"), figA,
       width=9, height=6.5)
message("Fig5a final saved")
