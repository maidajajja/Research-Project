suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

df <- read.csv("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/gene_location/convergent_strain_summary.csv",
               stringsAsFactors=FALSE)

# Parse plasmid columns from replicons string
df$IncFIB <- as.integer(grepl("IncFIB", df$replicons))
df$IncFIA <- as.integer(grepl("IncFIA", df$replicons))
df$IncFII <- as.integer(grepl("IncFII", df$replicons))

# Sort ST11 by: IncHI1B desc → total virulence loci desc → Res desc
vir_cols <- c("Aerobactin_bin","rmpA2_bin","Yersiniabactin_bin",
              "Salmochelin_bin","RmpADC_bin","Colibactin_bin")
df$n_vir_loci <- rowSums(df[, vir_cols])

df$ST_order <- factor(df$ST, levels=c("ST11","ST23","ST6086"))
df <- df %>%
  arrange(ST_order, desc(n_vir_loci), desc(IncHI1B), desc(Res))

st_counts <- table(df$ST_order)
st_labels <- paste0(names(st_counts), " (n=", st_counts, ")")
names(st_labels) <- names(st_counts)

# Virulence loci matrix — reordered per feedback
gene_cols   <- c("Aerobactin_bin","Yersiniabactin_bin","Salmochelin_bin",
                 "rmpA2_bin","RmpADC_bin","Colibactin_bin")
gene_labels <- c("Aerobactin","Yersiniabactin","Salmochelin",
                 "rmpA2","RmpADC","Colibactin")
mat <- as.matrix(df[, gene_cols])
colnames(mat) <- gene_labels
rownames(mat) <- df$strain

st_cols  <- c("ST11"="#225522","ST23"="#774411","ST6086"="#BBBBBB")
bin_cols <- c("0"="grey92","1"="#6699CC")

# Row annotation — reordered: ST, Resistance, Virulence, plasmids
row_ann <- rowAnnotation(
  ST          = df$ST,
  Virulence   = df$Vir,
  Resistance  = df$Res,
  IncHI1B     = as.character(df$IncHI1B),
  IncFIB      = as.character(df$IncFIB),
  IncFIA      = as.character(df$IncFIA),
  IncFII      = as.character(df$IncFII),
  ESBL        = as.character(df$ESBL),
  Carbapenemase = as.character(df$Carbapenemase),
  col = list(
    ST         = st_cols,
    Virulence  = colorRamp2(c(4,5), c("#EEBB88","#CC6677")),
    Resistance = colorRamp2(c(2,3), c("#88BBDD","#332288")),
    IncHI1B    = bin_cols,
    IncFIB     = bin_cols,
    IncFIA     = bin_cols,
    IncFII     = bin_cols),
  annotation_label     = c("ST","Virulence","Resistance","IncHI1B","IncFIB","IncFIA","IncFII","ESBL","Carbapenemase"),
  annotation_name_gp   = gpar(fontsize=16, fontface="bold"),
  annotation_name_rot  = 90,
  annotation_name_side = "top",
  annotation_width     = unit(c(0.6,0.4,0.4,0.5,0.5,0.5,0.5,0.5,0.5),"cm"),
  gap = unit(c(1,1,5,1,1,1,5,1,1),"mm"),
  show_legend = c(FALSE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE),
  annotation_legend_param = list(
    ST = list(
      title="ST", labels=c("ST11","ST23","ST6086"), at=c("ST11","ST23","ST6086"),
      title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16)),
    Virulence = list(
      title="Virulence score",
      title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16)),
    Resistance = list(
      title="Resistance score",
      title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16)),
    IncHI1B = list(
      title="Replicon present", labels=c("Absent","Present"), at=c("0","1"),
      title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16)),
    ESBL = list(
      title="ESBL", labels=c("Absent","Present"), at=c("0","1"),
      title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16)),

    Carbapenemase = list(
      title="Carbapenemase", labels=c("Absent","Present"), at=c("0","1"),
      title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16))))

# Column labels block annotation
col_block <- HeatmapAnnotation(
  `Hypervirulence loci` = anno_block(
    gp = gpar(fill="#E8EEF4", col="grey60", lwd=0.8),
    labels = "Hypervirulence loci",
    labels_gp = gpar(fontsize=16, fontface="bold", col="grey20")),
  show_annotation_name = FALSE)


ht <- Heatmap(mat,
  name = "Virulence locus",
  col  = c("0"="grey96","1"="#4477AA"),
  left_annotation  = row_ann,
  top_annotation   = col_block,
  show_row_names   = FALSE,
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_gp  = gpar(fontsize=16, fontface="bold.italic"),
  column_names_rot = 35,
  column_names_side= "bottom",
  row_split  = df$ST_order,
  row_title  = " ",
  row_gap    = unit(4,"mm"),
  height     = unit(nrow(mat)*0.55,"cm"),
  row_title_gp  = gpar(fontsize=16, fontface="bold"),
  row_title_rot = 0,
  column_title  = NULL,
  rect_gp = gpar(col="white", lwd=1.5),
  heatmap_legend_param = list(
    title="Virulence locus",
    labels=c("Absent","Present"), at=c(0,1),
    legend_gp=gpar(fill=c("grey96","#4477AA")),
    title_gp=gpar(fontsize=16,fontface="bold"),
    labels_gp=gpar(fontsize=16)),
  border=TRUE,
  border_gp=gpar(col="grey85", lwd=0.3))

png(file.path(OUT,"Fig8_convergent_plasmid_heatmap_v4.png"),
    width=15, height=10, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(20,5,30,10),"mm"))

dev.off()

pdf(file.path(OUT,"Fig8_convergent_plasmid_heatmap_v4.pdf"), width=15, height=10)
draw(ht, merge_legend=TRUE, padding=unit(c(20,5,30,10),"mm"))
# Section headers above row annotation groups

dev.off()

message("Fig8 v4 saved")
