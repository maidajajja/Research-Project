suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

df <- read.csv("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/gene_location/convergent_strain_summary.csv",
               stringsAsFactors=FALSE)

# Order: ST23 first, then ST11 by IncHI1B status, then ST6086
df$ST_order <- factor(df$ST, levels=c("ST23","ST11","ST6086"))
df <- df[order(df$ST_order, -df$IncHI1B, -df$Vir),]

# Build gene matrix
gene_cols <- c("Aerobactin_bin","Salmochelin_bin","RmpADC_bin","rmpA2_bin","Yersiniabactin_bin","Colibactin_bin")
gene_labels <- c("Aerobactin","Salmochelin","RmpADC","rmpA2","Yersiniabactin","Colibactin")
mat <- as.matrix(df[, gene_cols])
colnames(mat) <- gene_labels
rownames(mat) <- df$strain

# Colour palettes — match existing figures
st_cols <- c("ST23"="#E69F00","ST11"="#009E73","ST6086"="#4D4D4D")
inchi_cols <- c("0"="grey88","1"="#7B3294")

row_ann <- rowAnnotation(
  ST = df$ST,
  `Vir Score` = df$Vir,
  `Res Score` = df$Res,
  `IncHI1B` = as.character(df$IncHI1B),
  col = list(
    ST = st_cols,
    `Vir Score` = colorRamp2(c(4,5), c("#FE9929","#CC4C02")),
    `Res Score` = colorRamp2(c(2,3), c("#66C2A4","#00441B")),
    `IncHI1B` = inchi_cols
  ),
  annotation_label = c("ST","Vir Score","Res Score","IncHI1B"),
  annotation_name_gp = gpar(fontsize=11, fontface="bold"),
  annotation_name_rot = 90,
  annotation_name_side = "top",
  annotation_width = unit(c(0.6,0.4,0.4,0.5),"cm"),
  gap = unit(2,"mm"),
  show_legend = c(TRUE,TRUE,TRUE,TRUE),
  annotation_legend_param = list(
    ST = list(title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=10)),
    `Vir Score` = list(title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=10)),
    `Res Score` = list(title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=10)),
    `IncHI1B` = list(
      title="IncHI1B plasmid",
      labels=c("Absent","Present"),
      at=c("0","1"),
      title_gp=gpar(fontsize=11,fontface="bold"),
      labels_gp=gpar(fontsize=10)
    )
  )
)

ht <- Heatmap(mat,
  name = "Virulence gene",
  col = c("0"="#F0F0F0","1"="#08519C"),
  left_annotation = row_ann,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize=13, fontface="italic"),
  column_names_rot = 35,
  column_names_side = "bottom",
  row_split = df$ST_order,
  row_gap = unit(4,"mm"),
  row_title_gp = gpar(fontsize=13, fontface="bold"),
  row_title_rot = 0,
  column_title = NULL,
  rect_gp = gpar(col="white", lwd=1),
  heatmap_legend_param = list(
    title="Virulence gene",
    labels=c("Absent","Present"),
    at=c(0,1),
    legend_gp=gpar(fill=c("#F0F0F0","#08519C")),
    title_gp=gpar(fontsize=11,fontface="bold"),
    labels_gp=gpar(fontsize=10)
  )
)

png(file.path(OUT,"Fig8_convergent_plasmid_heatmap.png"),
    width=10, height=8, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,30,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig8_convergent_plasmid_heatmap.pdf"), width=10, height=8)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,30,5),"mm"))
dev.off()

message("Fig8 saved")
