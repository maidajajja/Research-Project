suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Source  <- meta[["Isolation Source"]]
meta$Health  <- meta[["Host Health"]]
meta$Source[is.na(meta$Source)|meta$Source==""] <- "Unknown"
meta$Health[is.na(meta$Health)|meta$Health==""] <- "Unknown"
meta <- meta[meta$Sample!="", c("Sample","Source","Health")]

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

st_col  <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]

kleb$ST  <- gsub("-.*","", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))

# Exact virulence determinant columns
vir_cols <- c(
  "klebsiella__ybst__Yersiniabactin",
  "klebsiella__cbst__Colibactin",
  "klebsiella__abst__Aerobactin",
  "klebsiella__smst__Salmochelin",
  "klebsiella__rmst__RmpADC",
  "klebsiella__rmpa2__rmpA2"
)
vir_cols <- vir_cols[vir_cols %in% colnames(kleb)]
message("Using virulence cols: ", paste(vir_cols, collapse=", "))

st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST,
                     Vir=kleb$Vir, Res=kleb$Res,
                     stringsAsFactors=FALSE)

# Build binary matrix
vir_df <- kleb[, c("strain", vir_cols)]
mat <- matrix(0, nrow=nrow(vir_df), ncol=length(vir_cols),
              dimnames=list(vir_df$strain, vir_cols))
for (col in vir_cols) {
  mat[, col] <- ifelse(is.na(vir_df[[col]]) |
                       vir_df[[col]] == "-" |
                       vir_df[[col]] == "0" |
                       vir_df[[col]] == "", 0, 1)
}

# Clean column names for display
colnames(mat) <- c("Yersiniabactin","Colibactin","Aerobactin",
                   "Salmochelin","RmpADC","rmpA2")[1:ncol(mat)]

ann <- data.frame(Sample=rownames(mat), stringsAsFactors=FALSE) %>%
  left_join(st_map, by="Sample") %>%
  left_join(meta, by="Sample")
ann$ST[is.na(ann$ST)]         <- "Unknown"
ann$Source[is.na(ann$Source)] <- "Unknown"
ann$Health[is.na(ann$Health)] <- "Unknown"
ann$Vir[is.na(ann$Vir)]       <- 0
ann$Res[is.na(ann$Res)]       <- 0

top_sts <- names(sort(table(ann$ST[ann$ST!="Unknown"]), decreasing=TRUE))[1:7]
st_pal  <- brewer.pal(7,"Set1")
st_cols <- setNames(c(st_pal,"grey80","grey90"), c(top_sts,"Other","Unknown"))
ann$ST_group <- ifelse(ann$ST %in% top_sts, ann$ST,
                       ifelse(ann$ST=="Unknown","Unknown","Other"))
ann$ST_split <- factor(ann$ST_group, levels=c(top_sts[1:5],"Other","Unknown"))

top_sources <- names(sort(table(ann$Source), decreasing=TRUE))[1:5]
ann$Src_group <- ifelse(ann$Source %in% top_sources, ann$Source, "Other")
src_cols <- setNames(brewer.pal(max(length(unique(ann$Src_group)),3),"Set2")[1:length(unique(ann$Src_group))],
                     unique(ann$Src_group))

health_vals <- unique(ann$Health)
health_cols <- setNames(colorRampPalette(brewer.pal(8,"Accent"))(length(health_vals)),
                        health_vals)

row_ann <- rowAnnotation(
  ST           = ann$ST_group,
  `Vir. Score` = ann$Vir,
  Source       = ann$Src_group,
  col = list(
    ST           = st_cols,
    `Vir. Score` = colorRamp2(c(0,2,5), c("#F7FBFF","#6BAED6","#08306B")),
    Source       = src_cols
  ),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_width   = unit(c(0.8,0.4,0.6),"cm"),
  gap = unit(1,"mm")
)

ht <- Heatmap(mat,
  name  = "Virulence",
  col   = c("0"="#F0F0F0","1"="#08519C"),
  left_annotation  = row_ann,
  show_row_names   = FALSE,
  column_names_gp  = gpar(fontsize=10, fontface="italic"),
  column_names_rot = 45,
  row_split        = ann$ST_split,
  row_gap          = unit(3,"mm"),
  row_title_gp     = gpar(fontsize=9, fontface="bold"),
  row_title_rot    = 0,
  column_title     = expression(italic("K. pneumoniae")~"Virulence Determinants (n = 234)"),
  column_title_gp  = gpar(fontsize=12, fontface="bold"),
  clustering_distance_rows    = "binary",
  clustering_distance_columns = "binary",
  clustering_method_rows      = "ward.D2",
  clustering_method_columns   = "ward.D2",
  rect_gp = gpar(col="white", lwd=0.5),
  heatmap_legend_param = list(
    title="Present", labels=c("No","Yes"),
    at=c(0,1), legend_height=unit(2,"cm")
  )
)

png(file.path(OUT,"Fig3_virulence_heatmap.png"),
    width=9, height=14, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig3_virulence_heatmap.pdf"), width=9, height=14)
draw(ht, merge_legend=TRUE)
dev.off()

message("Fig3 saved")
