suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load Rtab
rtab <- read.table("/scratch/users/k22017808/KP_Research_Project/06_Pangenome_final229/gene_presence_absence.Rtab",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
rownames(rtab) <- rtab[,1]
rtab <- rtab[,-1]
mat <- as.matrix(rtab)  # genes x samples

n_samples <- ncol(mat)
n_genes <- nrow(mat)
freq <- rowSums(mat)

# Categorise genes
gene_cat <- ifelse(freq >= 0.99*n_samples, "Core",
            ifelse(freq >= 0.95*n_samples, "Soft core", "Shell"))

# Keep core + soft core + top variable shell genes
core_idx <- which(gene_cat %in% c("Core","Soft core"))
shell_idx <- which(gene_cat == "Shell")
shell_var <- apply(mat[shell_idx,], 1, var)
shell_top_idx <- shell_idx[order(shell_var, decreasing=TRUE)[1:300]]
keep_idx <- c(core_idx, shell_top_idx)

mat_filt <- mat[keep_idx,]
gene_cat_filt <- gene_cat[keep_idx]

# Order genes: core first, soft core, shell
gene_order <- order(factor(gene_cat_filt, levels=c("Core","Soft core","Shell")),
                    -rowSums(mat_filt))
mat_filt <- mat_filt[gene_order,]
gene_cat_filt <- gene_cat_filt[gene_order]

message(sprintf("Filtered matrix: %d genes x %d samples", nrow(mat_filt), ncol(mat_filt)))

# Load ST metadata
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample <- as.character(meta[["Genome ID"]])
meta$ST_meta <- sub(".*\\.","", meta[["MLST"]])

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
kleb$ST <- gsub("^ST","", kleb$ST)
st_map <- setNames(kleb$ST, kleb$strain)

# get ST per sample
sample_ids <- colnames(mat_filt)
sample_st <- sapply(sample_ids, function(s){
  st <- st_map[s]
  if(is.na(st)||st=="") st <- meta$ST_meta[meta$Sample==s][1]
  if(is.na(st)||st=="") st <- "Unknown"
  as.character(st)
})

top_sts <- names(sort(table(sample_st[sample_st!="Unknown"]),
                      decreasing=TRUE))[1:8]
sample_st_group <- ifelse(sample_st %in% top_sts,
                          paste0("ST",sample_st), "Other")

# order samples by ST then by total gene count
sample_order <- order(factor(sample_st_group,
                             levels=c("Other",rev(paste0("ST",top_sts)))),
                      -colSums(mat_filt))
mat_ord <- mat_filt[, sample_order]
sample_st_ord <- sample_st_group[sample_order]

# TRANSPOSE: samples as rows, genes as columns
mat_t <- t(mat_ord)

# CBF colours
cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#4D4D4D","grey75")
st_levels <- c(paste0("ST",top_sts),"Other")
st_cols <- setNames(cbf_pal[1:length(st_levels)], st_levels)

cat_cols <- c("Core"="#08519C","Soft core"="#6BAED6","Shell"="#E69F00")

# Row annotation (isolates = rows)
row_ann <- rowAnnotation(
  ST = sample_st_ord,
  col = list(ST = st_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  simple_anno_size = unit(5,"mm"),
  show_legend = TRUE
)

# Column annotation (genes = columns)
col_ann <- HeatmapAnnotation(
  Category = gene_cat_filt[gene_order],
  col = list(Category = cat_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  simple_anno_size = unit(4,"mm"),
  show_legend = TRUE
)

ht <- Heatmap(mat_t,
  name = "Gene",
  col = c("0"="grey97","1"="#08519C"),
  left_annotation = row_ann,
  top_annotation = col_ann,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_split = factor(sample_st_ord, levels=c(paste0("ST",top_sts),"Other")),
  row_title_gp = gpar(fontsize=8, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(1,"mm"),
  column_split = factor(gene_cat_filt[gene_order],
                        levels=c("Core","Soft core","Shell")),
  column_title_gp = gpar(fontsize=9, fontface="bold"),
  column_gap = unit(2,"mm"),
  heatmap_legend_param = list(
    title = "Gene",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey97","#08519C")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)
  ),
  border = TRUE,
  use_raster = FALSE
)

png(file.path(OUT,"Fig6c_pangenome_heatmap.png"),
    width=16, height=10, units="in", res=600, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig6c_pangenome_heatmap.pdf"), width=16, height=10)
draw(ht, merge_legend=TRUE)
dev.off()

message("Fig6c v2 saved at 600 DPI")
