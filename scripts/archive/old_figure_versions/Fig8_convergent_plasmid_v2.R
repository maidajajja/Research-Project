suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

df <- read.csv("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/gene_location/convergent_strain_summary.csv",
               stringsAsFactors=FALSE)

cat("Convergent isolates by ST:\n")
print(table(df$ST))
cat("Columns:\n"); print(colnames(df))

# Order: ST11 first (most convergent), then ST23, then ST6086
df$ST_order <- factor(df$ST, levels=c("ST11","ST23","ST6086"))
df <- df[order(df$ST_order, -df$IncHI1B, -df$Vir), ]

# ST labels with n
st_counts <- table(df$ST_order)
st_labels <- paste0(names(st_counts), " (n=", st_counts, ")")
names(st_labels) <- names(st_counts)

# Gene matrix
gene_cols <- c("Aerobactin_bin","Salmochelin_bin","RmpADC_bin",
               "rmpA2_bin","Yersiniabactin_bin","Colibactin_bin")
gene_labels <- c("Aerobactin","Salmochelin","RmpADC","rmpA2","Yersiniabactin","Colibactin")
mat <- as.matrix(df[, gene_cols])
colnames(mat) <- gene_labels
rownames(mat) <- df$strain

# Paul Tol muted colours — consistent with all other figures
st_cols <- c(
  "ST11"  = "#225522",  # dark forest green
  "ST23"  = "#774411",  # dark brown
  "ST6086"= "#BBBBBB"   # light grey
)

inchi_cols <- c("0"="grey92", "1"="#6699CC")  # muted blue not purple

row_ann <- rowAnnotation(
  ST = df$ST,
  `Vir Score` = df$Vir,
  `Res Score` = df$Res,
  `IncHI1B` = as.character(df$IncHI1B),
  col = list(
    ST = st_cols,
    `Vir Score` = colorRamp2(c(4,5), c("#EEBB88","#CC6677")),
    `Res Score` = colorRamp2(c(2,3), c("#88BBDD","#332288")),
    `IncHI1B` = inchi_cols),
  annotation_label = c("ST","Vir Score","Res Score","IncHI1B plasmid"),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_name_rot = 90,
  annotation_name_side = "top",
  annotation_width = unit(c(0.6,0.4,0.4,0.5),"cm"),
  gap = unit(2,"mm"),
  annotation_legend_param = list(
    ST = list(
      title="ST",
      labels=c("ST11","ST23","ST6086"),
      at=c("ST11","ST23","ST6086"),
      title_gp=gpar(fontsize=9,fontface="bold"),
      labels_gp=gpar(fontsize=8)),
    `Vir Score` = list(
      title="Virulence score",
      title_gp=gpar(fontsize=9,fontface="bold"),
      labels_gp=gpar(fontsize=8)),
    `Res Score` = list(
      title="Resistance score",
      title_gp=gpar(fontsize=9,fontface="bold"),
      labels_gp=gpar(fontsize=8)),
    `IncHI1B` = list(
      title="IncHI1B plasmid",
      labels=c("Absent","Present"),
      at=c("0","1"),
      title_gp=gpar(fontsize=9,fontface="bold"),
      labels_gp=gpar(fontsize=8))))

ht <- Heatmap(mat,
  name = "Virulence locus",
  col = c("0"="grey96","1"="#4477AA"),
  left_annotation = row_ann,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize=11, fontface="bold.italic"),
  column_names_rot = 35,
  column_names_side = "bottom",
  row_split = df$ST_order,
  row_title = st_labels[levels(df$ST_order)],
  row_gap = unit(4,"mm"),
  row_title_gp = gpar(fontsize=11, fontface="bold"),
  row_title_rot = 0,
  column_title = NULL,
  rect_gp = gpar(col="white", lwd=1.5),
  heatmap_legend_param = list(
    title="Virulence locus",
    labels=c("Absent","Present"),
    at=c(0,1),
    legend_gp=gpar(fill=c("grey96","#4477AA")),
    title_gp=gpar(fontsize=9,fontface="bold"),
    labels_gp=gpar(fontsize=8)),
  border=TRUE,
  border_gp=gpar(col="grey85", lwd=0.3))

png(file.path(OUT,"Fig8_convergent_plasmid_heatmap_v2.png"),
    width=10, height=8, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,30,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig8_convergent_plasmid_heatmap_v2.pdf"), width=10, height=8)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,30,5),"mm"))
dev.off()

message("Fig8 v2 saved")
