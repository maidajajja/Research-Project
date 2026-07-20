suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load data
df <- read.csv("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/gene_location/convergent_strain_summary.csv",
               stringsAsFactors=FALSE)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

# Parse plasmid columns
df$IncFIB <- as.integer(grepl("IncFIB", df$replicons))
df$IncFIA <- as.integer(grepl("IncFIA", df$replicons))
df$IncFII <- as.integer(grepl("IncFII", df$replicons))

# Add AMR class columns from Kleborate
amr_map <- list(
  AGly = "klebsiella_pneumo_complex__amr__AGly_acquired",
  Flq  = "klebsiella_pneumo_complex__amr__Flq_acquired",
  Sul  = "klebsiella_pneumo_complex__amr__Sul_acquired",
  Tet  = "klebsiella_pneumo_complex__amr__Tet_acquired",
  Phe  = "klebsiella_pneumo_complex__amr__Phe_acquired",
  Tmt  = "klebsiella_pneumo_complex__amr__Tmt_acquired",
  MLS  = "klebsiella_pneumo_complex__amr__MLS_acquired",
  Rif  = "klebsiella_pneumo_complex__amr__Rif_acquired"
)
for (nm in names(amr_map)) {
  col <- amr_map[[nm]]
  df[[nm]] <- sapply(df$strain, function(s) {
    val <- kleb[kleb$strain == s, col]
    if (length(val)==0 || val=="-") 0L else 1L
  })
}

# Sort: ST group → n_vir_loci desc → IncHI1B desc → Res desc
vir_cols <- c("Aerobactin_bin","Yersiniabactin_bin","Salmochelin_bin",
              "rmpA2_bin","RmpADC_bin","Colibactin_bin")
df$n_vir_loci <- rowSums(df[, vir_cols])
df$ST_order <- factor(df$ST, levels=c("ST11","ST23","ST6086"))
df <- df %>% arrange(ST_order, desc(n_vir_loci), desc(IncHI1B), desc(Res))

st_counts <- table(df$ST_order)
st_labels <- paste0(names(st_counts), " (n=", st_counts, ")")
names(st_labels) <- names(st_counts)

# Colour schemes
bin_cols    <- c("0"="grey92", "1"="#6699CC")
esbl_cols   <- c("0"="grey92", "1"="#44AA99")
carb_cols   <- c("0"="grey92", "1"="#AA3377")
amr_cols_cl <- c("0"="grey92", "1"="#555555")
st_cols     <- c("ST11"="#225522","ST23"="#774411","ST6086"="#BBBBBB")

# Left annotation: metadata + plasmids + ESBL + Carb
row_ann <- rowAnnotation(
  ST      = df$ST,
  IncHI1B = as.character(df$IncHI1B),
  IncFIB  = as.character(df$IncFIB),
  IncFIA  = as.character(df$IncFIA),
  IncFII  = as.character(df$IncFII),
  col = list(
    ST      = st_cols,
    IncHI1B = bin_cols,
    IncFIB  = bin_cols,
    IncFIA  = bin_cols,
    IncFII  = bin_cols),
  annotation_label     = c("ST","IncHI1B","IncFIB","IncFIA","IncFII"),
  annotation_name_gp   = gpar(fontsize=13, fontface="bold"),
  annotation_name_rot  = 90,
  annotation_name_side = "top",
  annotation_width     = unit(c(0.6,0.5,0.5,0.5,0.5),"cm"),
  gap = unit(c(1,5,1,1,1),"mm"),
  show_legend = c(TRUE,TRUE,FALSE,FALSE,FALSE),
  annotation_legend_param = list(
    ST = list(title="ST", labels=c("ST11","ST23","ST6086"), at=c("ST11","ST23","ST6086"),
      title_gp=gpar(fontsize=13,fontface="bold"), labels_gp=gpar(fontsize=12)),
    IncHI1B = list(title="Replicon present", labels=c("Absent","Present"), at=c("0","1"),
      title_gp=gpar(fontsize=13,fontface="bold"), labels_gp=gpar(fontsize=12))))

# Block 1: Virulence loci matrix
vir_mat <- as.matrix(df[, vir_cols])
colnames(vir_mat) <- c("Aerobactin","Yersiniabactin","Salmochelin","rmpA2","RmpADC","Colibactin")
rownames(vir_mat) <- df$strain

# Block 2: AMR classes matrix
amr_class_cols <- c("ESBL","Carbapenemase","AGly","Flq","Sul","Tet","Phe","Tmt","MLS","Rif")
amr_labels     <- c("ESBL","Carbapenemase","Aminoglycoside","Fluoroquinolone",
                    "Sulphonamide","Tetracycline","Phenicol","Trimethoprim","Macrolide","Rifampicin")
amr_mat <- as.matrix(df[, amr_class_cols])
amr_mat <- apply(amr_mat, 2, as.integer)
colnames(amr_mat) <- amr_labels
rownames(amr_mat) <- df$strain

# Top annotations for each block
top_vir <- HeatmapAnnotation(
  ` ` = anno_block(gp=gpar(fill="#E8EEF4", col="grey60", lwd=0.8),
    labels="Hypervirulence loci",
    labels_gp=gpar(fontsize=13, fontface="bold", col="grey20")),
  show_annotation_name=FALSE)

top_amr <- HeatmapAnnotation(
  ` ` = anno_block(gp=gpar(fill="#F5E8EA", col="grey60", lwd=0.8),
    labels="AMR determinants",
    labels_gp=gpar(fontsize=13, fontface="bold", col="grey20")),
  show_annotation_name=FALSE)

ht_vir <- Heatmap(vir_mat,
  name = "Virulence locus",
  col  = c("0"="grey96","1"="#4477AA"),
  top_annotation   = top_vir,
  left_annotation  = row_ann,
  show_row_names   = FALSE,
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_gp  = gpar(fontsize=13, fontface="bold.italic"),
  column_names_rot = 35,
  column_names_side= "bottom",
  row_split  = df$ST_order,
  row_title  = " ",
  row_gap    = unit(8,"mm"),
  row_title_gp = gpar(fontsize=13, fontface="bold"),
  column_title = NULL,
  rect_gp = gpar(col="white", lwd=1.5),
  height = unit(nrow(vir_mat)*0.55,"cm"),
  heatmap_legend_param = list(
    title="Virulence locus", labels=c("Absent","Present"), at=c(0,1),
    legend_gp=gpar(fill=c("grey96","#4477AA")),
    title_gp=gpar(fontsize=13,fontface="bold"), labels_gp=gpar(fontsize=12)),
  border=TRUE, border_gp=gpar(col="grey85", lwd=0.3))

ht_amr <- Heatmap(amr_mat,
  name = "AMR class",
  col  = c("0"="grey96","1"="#555555"),
  top_annotation   = top_amr,
  show_row_names   = FALSE,
  show_row_dend    = FALSE,
  show_column_dend = FALSE,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_names_gp  = gpar(fontsize=13, fontface="bold.italic"),
  column_names_rot = 35,
  column_names_side= "bottom",
  row_split  = df$ST_order,
  row_title  = " ",
  row_gap    = unit(8,"mm"),
  column_title = NULL,
  rect_gp = gpar(col="white", lwd=1.5),
  height = unit(nrow(amr_mat)*0.55,"cm"),
  heatmap_legend_param = list(
    title="AMR class", labels=c("Absent","Present"), at=c(0,1),
    legend_gp=gpar(fill=c("grey96","#555555")),
    title_gp=gpar(fontsize=13,fontface="bold"), labels_gp=gpar(fontsize=12)),
  border=TRUE, border_gp=gpar(col="grey85", lwd=0.3))

png(file.path(OUT,"Fig8b_convergent_integrated_v1.png"),
    width=17, height=10, units="in", res=300, bg="white")
draw(ht_vir + ht_amr, merge_legend=TRUE,
     padding=unit(c(20,5,30,10),"mm"))
dev.off()

pdf(file.path(OUT,"Fig8b_convergent_integrated_v1.pdf"), width=17, height=10)
draw(ht_vir + ht_amr, merge_legend=TRUE,
     padding=unit(c(20,5,30,10),"mm"))
dev.off()

message("Fig8b v1 saved")
