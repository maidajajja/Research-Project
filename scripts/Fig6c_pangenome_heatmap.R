suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

rtab <- read.table("/scratch/users/k22017808/KP_Research_Project/06_Pangenome/gene_presence_absence.Rtab",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
rownames(rtab) <- rtab[,1]
rtab <- rtab[,-1]
mat <- as.matrix(rtab)

n <- ncol(mat)
freq <- rowSums(mat)

core_soft <- mat[freq >= 0.95*n, ]
shell <- mat[freq >= 0.15*n & freq < 0.95*n, ]
shell_var <- apply(shell, 1, var)
shell_top <- shell[order(shell_var, decreasing=TRUE)[1:200], ]
mat_filt <- rbind(core_soft, shell_top)
message(sprintf("Matrix: %d genes x %d isolates", nrow(mat_filt), ncol(mat_filt)))

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample <- as.character(meta[["Genome ID"]])

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
st_map <- setNames(kleb$ST, kleb$strain)

sample_ids <- colnames(mat_filt)
sample_st <- sapply(sample_ids, function(s){
  st <- st_map[s]
  if(is.na(st)) st <- "Unknown"
  st
})

top_sts <- names(sort(table(sample_st[sample_st!="Unknown"]), decreasing=TRUE))[1:8]
sample_st_group <- ifelse(sample_st %in% top_sts, sample_st, "Other")
sample_order <- order(factor(sample_st_group, levels=c(top_sts,"Other")))
mat_ord <- mat_filt[, sample_order]
sample_st_ord <- sample_st_group[sample_order]

gene_freq <- rowSums(mat_filt)
gene_cat <- ifelse(gene_freq >= 0.99*n, "Core",
            ifelse(gene_freq >= 0.95*n, "Soft core", "Shell"))

wong_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#4D4D4D","grey75")
st_levels <- c(top_sts,"Other")
st_cols <- setNames(wong_pal[1:length(st_levels)], st_levels)

cat_cols <- c("Core"="#08306B","Soft core"="#4292C6","Shell"="#FC8D59")

col_ann <- HeatmapAnnotation(
  ST = sample_st_ord,
  col = list(ST = st_cols),
  annotation_name_gp = gpar(fontsize=12, fontface="bold"),
  annotation_label = "ST",
  simple_anno_size = unit(5,"mm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize=12, fontface="bold"),
    labels_gp = gpar(fontsize=11)
  )
)

row_ann <- rowAnnotation(
  Category = gene_cat,
  col = list(Category = cat_cols),
  annotation_name_gp = gpar(fontsize=12, fontface="bold"),
  annotation_label = "Category",
  simple_anno_size = unit(5,"mm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize=12, fontface="bold"),
    labels_gp = gpar(fontsize=11)
  )
)

ht <- Heatmap(mat_ord,
  name = "Gene presence",
  col = c("0"="grey95","1"="#08519C"),
  top_annotation = col_ann,
  right_annotation = row_ann,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = sample_st_ord,
  column_title_gp = gpar(fontsize=11, fontface="bold"),
  column_title_rot = 90,
  column_gap = unit(1,"mm"),
  heatmap_legend_param = list(
    title = "Gene",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey95","#08519C")),
    title_gp = gpar(fontsize=12, fontface="bold"),
    labels_gp = gpar(fontsize=11)
  ),
  border = TRUE,
  use_raster = FALSE
)

png(file.path(OUT,"Fig6c_pangenome_heatmap.png"),
    width=14, height=10, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig6c_pangenome_heatmap.pdf"), width=14, height=10)
draw(ht, merge_legend=TRUE)
dev.off()

message("Fig6c saved")
