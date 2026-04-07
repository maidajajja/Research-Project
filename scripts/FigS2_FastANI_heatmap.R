suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

ani <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI/fastani_results.txt",
                  sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(ani) <- c("Query","Reference","ANI","Fragments","Total")

clean_name <- function(x) {
  x <- basename(x)
  x <- gsub(".fna$","", x)
  x <- gsub("_genomic$","", x)
  x
}
ani$Query     <- clean_name(ani$Query)
ani$Reference <- clean_name(ani$Reference)

all_samples <- unique(c(ani$Query, ani$Reference))
n <- length(all_samples)
message("Samples: ", n)

mat <- matrix(NA, nrow=n, ncol=n, dimnames=list(all_samples, all_samples))
diag(mat) <- 100
for (i in 1:nrow(ani)) {
  mat[ani$Query[i], ani$Reference[i]] <- ani$ANI[i]
  mat[ani$Reference[i], ani$Query[i]] <- ani$ANI[i]
}
mat[is.na(mat)] <- min(ani$ANI, na.rm=TRUE)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample <- as.character(meta[["Genome ID"]])
meta$Country <- meta[["Isolation Country"]]
meta$Country[is.na(meta$Country)|meta$Country==""] <- "Unknown"
meta <- meta[, c("Sample","Country")]

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, stringsAsFactors=FALSE)

ann_df <- data.frame(Sample=rownames(mat), stringsAsFactors=FALSE) %>%
  left_join(st_map, by="Sample") %>%
  left_join(meta, by="Sample")
ann_df$ST[is.na(ann_df$ST)] <- "Unknown"
ann_df$Country[is.na(ann_df$Country)] <- "Unknown"

top_sts <- names(sort(table(ann_df$ST[ann_df$ST!="Unknown"]), decreasing=TRUE))[1:7]
ann_df$ST_group <- ifelse(ann_df$ST %in% top_sts, ann_df$ST, "Other")
top_countries <- names(sort(table(ann_df$Country), decreasing=TRUE))[1:6]
ann_df$Cty_group <- ifelse(ann_df$Country %in% top_countries, ann_df$Country, "Other")

ann_row <- data.frame(ST=ann_df$ST_group, Country=ann_df$Cty_group, row.names=ann_df$Sample)
st_pal <- brewer.pal(7,"Set1")
st_cols <- setNames(c(st_pal,"grey80"), c(top_sts,"Other"))
cty_cols <- setNames(colorRampPalette(brewer.pal(7,"Set3"))(length(unique(ann_df$Cty_group))), unique(ann_df$Cty_group))
ann_colors <- list(ST=st_cols, Country=cty_cols)

ani_cols <- colorRampPalette(c("#3D0357","#4B0082","#0D3B8E","#008080","#00CED1","#FFFF00"))(100)

png(file.path(OUT,"FigS2_FastANI_heatmap.png"), width=14, height=12, units="in", res=300, bg="white")
pheatmap(mat, color=ani_cols, annotation_row=ann_row, annotation_col=ann_row,
         annotation_colors=ann_colors, show_rownames=FALSE, show_colnames=FALSE,
         clustering_method="ward.D2", border_color=NA,
         main="K. pneumoniae Pairwise Average Nucleotide Identity (n = 234)",
         fontsize=11, annotation_legend=TRUE)
dev.off()

pdf(file.path(OUT,"FigS2_FastANI_heatmap.pdf"), width=14, height=12)
pheatmap(mat, color=ani_cols, annotation_row=ann_row, annotation_col=ann_row,
         annotation_colors=ann_colors, show_rownames=FALSE, show_colnames=FALSE,
         clustering_method="ward.D2", border_color=NA,
         main="K. pneumoniae Pairwise Average Nucleotide Identity (n = 234)",
         fontsize=11, annotation_legend=TRUE)
dev.off()
message("FigS2 saved")
