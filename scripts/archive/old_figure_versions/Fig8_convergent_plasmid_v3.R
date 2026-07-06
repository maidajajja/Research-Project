suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

df <- read.csv("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/gene_location/convergent_strain_summary.csv",
               stringsAsFactors=FALSE)

# Parse plasmid binary columns from replicons string
df$IncFIB <- as.integer(grepl("IncFIB", df$replicons))
df$IncFIA <- as.integer(grepl("IncFIA", df$replicons))
df$IncFII <- as.integer(grepl("IncFII", df$replicons))

# Order rows
df$ST_order <- factor(df$ST, levels=c("ST11","ST23","ST6086"))
df <- df[order(df$ST_order, -df$IncHI1B, -df$Vir), ]

st_counts <- table(df$ST_order)
st_labels <- paste0(names(st_counts), " (n=", st_counts, ")")
names(st_labels) <- names(st_counts)

# Virulence loci matrix
gene_cols   <- c("Aerobactin_bin","Salmochelin_bin","RmpADC_bin","rmpA2_bin","Yersiniabactin_bin","Colibactin_bin")
gene_labels <- c("Aerobactin","Salmochelin","RmpADC","rmpA2","Yersiniabactin","Colibactin")
mat <- as.matrix(df[, gene_cols])
colnames(mat) <- gene_labels
rownames(mat) <- df$strain

st_cols <- c("ST11"="#225522","ST23"="#774411","ST6086"="#BBBBBB")
bin_cols <- c("0"="grey92","1"="#6699CC")

row_ann <- rowAnnotation(
  ST         = df$ST,
  `Vir Score`= df$Vir,
  `Res Score`= df$Res,
  IncHI1B    = as.character(df$IncHI1B),
  IncFIB     = as.character(df$IncFIB),
  IncFIA     = as.character(df$IncFIA),
  IncFII     = as.character(df$IncFII),
  col = list(
    ST          = st_cols,
    `Vir Score` = colorRamp2(c(4,5), c("#EEBB88","#CC6677")),
    `Res Score` = colorRamp2(c(2,3), c("#88BBDD","#332288")),
    IncHI1B     = bin_cols,
    IncFIB      = bin_cols,
    IncFIA      = bin_cols,
    IncFII      = bin_cols),
  annotation_label = c("ST","Vir Score","Res Score","IncHI1B","IncFIB","IncFIA","IncFII"),
  annotation_name_gp   = gpar(fontsize=12, fontface="bold"),
  annotation_name_rot  = 90,
  annotation_name_side = "top",
  annotation_width = unit(c(0.6,0.4,0.4,0.5,0.5,0.5,0.5),"cm"),
  gap = unit(2,"mm"),
  show_legend = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE),
  annotation_legend_param = list(
    ST = list(title="ST", labels=c("ST11","ST23","ST6086"), at=c("ST11","ST23","ST6086"),
              title_gp=gpar(fontsize=12,fontface="bold"), labels_gp=gpar(fontsize=11)),
    `Vir Score` = list(title="Virulence score",
              title_gp=gpar(fontsize=12,fontface="bold"), labels_gp=gpar(fontsize=11)),
    `Res Score` = list(title="Resistance score",
              title_gp=gpar(fontsize=12,fontface="bold"), labels_gp=gpar(fontsize=11)),
    IncHI1B = list(title="Plasmid replicon", labels=c("Absent","Present"), at=c("0","1"),
              title_gp=gpar(fontsize=12,fontface="bold"), labels_gp=gpar(fontsize=11))))

ht <- Heatmap(mat,
  name = "Virulence locus",
  col  = c("0"="grey96","1"="#4477AA"),
  left_annotation = row_ann,
  show_row_names   = FALSE,
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_gp  = gpar(fontsize=12, fontface="bold.italic"),
  column_names_rot = 35,
  column_names_side= "bottom",
  row_split  = df$ST_order,
  row_title  = st_labels[levels(df$ST_order)],
  row_gap    = unit(4,"mm"),
  row_title_gp  = gpar(fontsize=12, fontface="bold"),
  row_title_rot = 0,
  column_title  = NULL,
  rect_gp = gpar(col="white", lwd=1.5),
  heatmap_legend_param = list(
    title="Virulence locus", labels=c("Absent","Present"), at=c(0,1),
    legend_gp=gpar(fill=c("grey96","#4477AA")),
    title_gp=gpar(fontsize=12,fontface="bold"), labels_gp=gpar(fontsize=11)),
  border=TRUE, border_gp=gpar(col="grey85", lwd=0.3))

png(file.path(OUT,"Fig8_convergent_plasmid_heatmap_v3.png"),
    width=10, height=8, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,30,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig8_convergent_plasmid_heatmap_v3.pdf"), width=10, height=8)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,30,5),"mm"))
dev.off()

message("Fig8 v3 saved")
