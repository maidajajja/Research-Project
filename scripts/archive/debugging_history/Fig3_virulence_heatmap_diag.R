suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

# Build a Sample ID that actually matches Kleborate's strain column.
# genomes.csv only ever stores BV-BRC numeric IDs in "Genome ID", but
# Kleborate's strain names are a mix of BV-BRC IDs, GCA/GCF accessions,
# and SRA-based names depending on each genome's original source.
# Fall back to Assembly Accession (prefix match) then SRA Accession
# (substring match) when the raw Genome ID doesn't match directly.
meta$GenomeID_str <- as.character(meta[["Genome ID"]])
meta$Sample <- ifelse(meta$GenomeID_str %in% kleb$strain, meta$GenomeID_str, NA)
unmatched_idx <- which(is.na(meta$Sample))
for (i in unmatched_idx) {
  acc <- meta[["Assembly Accession"]][i]
  sra <- meta[["SRA Accession"]][i]
  if (!is.na(acc) && acc != "") {
    m <- kleb$strain[startsWith(kleb$strain, acc)]
    if (length(m) >= 1) meta$Sample[i] <- m[1]
  } else if (!is.na(sra) && sra != "") {
    m <- kleb$strain[grepl(sra, kleb$strain, fixed = TRUE)]
    if (length(m) >= 1) meta$Sample[i] <- m[1]
  }
}

meta$Source <- meta[["Isolation Source"]]
meta$Source[is.na(meta$Source)|meta$Source==""] <- "Unknown"
meta$Source <- dplyr::recode(meta$Source,
  "liver abscess" = "Abscess/Pus",
  "Liver abscess puncture fluid" = "Abscess/Pus",
  "liver abscess drainage fluid" = "Abscess/Pus",
  "Aspirate liver abscess" = "Abscess/Pus",
  "Drainage samples from liver abscess" = "Abscess/Pus",
  "Abscess drainage liquid" = "Abscess/Pus",
  "Abscess drainage" = "Abscess/Pus",
  "Pyogenic liver abscess" = "Abscess/Pus",
  "abscess hepatic" = "Abscess/Pus",
  "liver abscess of a 14-year-old boy" = "Abscess/Pus",
  "isolated from liver abcess in human, Buffalo, New York" = "Abscess/Pus",
  "Pus" = "Abscess/Pus",
  "pus" = "Abscess/Pus",
  "liver" = "Liver tissue/fluid",
  "Liver puncture" = "Liver tissue/fluid",
  "liver puncture fluid" = "Liver tissue/fluid",
  "Liver cyst" = "Liver tissue/fluid",
  "Liver drainage" = "Liver tissue/fluid",
  "Hepato/biliary" = "Liver tissue/fluid",
  "Bile" = "Liver tissue/fluid",
  "Bile fluid" = "Liver tissue/fluid",
  "Drainage fluid" = "Liver tissue/fluid",
  "Percutaneous drainage" = "Liver tissue/fluid",
  "abdominal drainage" = "Liver tissue/fluid",
  "abdominal fluid" = "Liver tissue/fluid",
  "Peritoneal fluid" = "Liver tissue/fluid",
  "A General Hospital in Southeast China" = "Liver tissue/fluid",
  "blood" = "Blood",
  "Blood" = "Blood",
  "blood, liver abscess" = "Blood",
  "feces" = "Gastrointestinal",
  "rectal swab" = "Gastrointestinal",
  "stool" = "Gastrointestinal",
  "Respiratory" = "Respiratory",
  "bronchoalveolar lavage" = "Respiratory",
  "patients" = "Unknown",
  "Physical" = "Unknown",
  "Other specimens" = "Unknown"
)
meta <- meta[!is.na(meta$Sample) & meta$Sample!="", c("Sample","Source")]

st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
kleb$ST <- gsub("-.*","", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, Vir=kleb$Vir, stringsAsFactors=FALSE)

vir_cols <- c(
  "klebsiella__ybst__Yersiniabactin",
  "klebsiella__cbst__Colibactin",
  "klebsiella__abst__Aerobactin",
  "klebsiella__smst__Salmochelin",
  "klebsiella__rmst__RmpADC",
  "klebsiella__rmpa2__rmpA2"
)
vir_cols <- vir_cols[vir_cols %in% colnames(kleb)]
vir_df <- kleb[, c("strain", vir_cols)]
mat <- matrix(0, nrow=nrow(vir_df), ncol=length(vir_cols),
              dimnames=list(vir_df$strain, vir_cols))
for (col in vir_cols) {
  mat[, col] <- ifelse(is.na(vir_df[[col]])|vir_df[[col]]=="-"|vir_df[[col]]=="0"|vir_df[[col]]=="", 0, 1)
}
colnames(mat) <- c("Yersiniabactin","Colibactin","Aerobactin","Salmochelin","RmpADC","rmpA2")[1:ncol(mat)]

ann <- data.frame(Sample=rownames(mat), stringsAsFactors=FALSE) %>%
  left_join(st_map, by="Sample") %>%
  left_join(meta, by="Sample")
ann$ST[is.na(ann$ST)] <- "Unknown"
ann$Source[is.na(ann$Source)] <- "Unknown"
ann$Vir[is.na(ann$Vir)] <- 0

cat("DIAGNOSTIC - ann$Source table right before plotting:\n")
print(table(ann$Source))
cat("\n")

top_sts <- names(sort(table(ann$ST[ann$ST!="Unknown"]), decreasing=TRUE))[1:7]
wong_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
st_cols <- setNames(c(wong_pal,"grey80"), c(top_sts,"Other"))
ann$ST_group <- ifelse(ann$ST %in% top_sts, ann$ST, "Other")
ann <- ann[!is.na(ann$ST_group), ]
mat <- mat[rownames(mat) %in% ann$Sample, ]
mat <- mat[order(rowSums(mat), decreasing=TRUE), ]
ann <- ann[match(rownames(mat), ann$Sample), ]

group_counts <- sort(table(ann$ST_group), decreasing=FALSE)
st_order <- names(group_counts)
ann$ST_split <- factor(ann$ST_group, levels=st_order)

src_cols <- c(
  "Abscess/Pus"        = "purple4",
  "Liver tissue/fluid" = "#C2A5CF",
  "Blood"              = "#D7191C",
  "Gastrointestinal"   = "#1B7837",
  "Respiratory"        = "#74ADD1",
  "Unknown"            = "grey85",
  "Other"              = "grey60"
)
ann$Src_group <- ifelse(ann$Source %in% names(src_cols), ann$Source, "Other")

row_ann <- rowAnnotation(
  ST          = ann$ST_group,
  `Vir Score` = ann$Vir,
  Source      = ann$Src_group,
  col = list(
    ST          = st_cols,
    `Vir Score` = colorRamp2(c(0,2,5), c("#FFF7BC","#FE9929","#CC4C02")),
    Source      = src_cols
  ),
  annotation_label     = c("ST","Vir Score","Source"),
  annotation_name_gp   = gpar(fontsize=11, fontface="bold"),
  annotation_name_rot  = 90,
  annotation_name_side = "top",
  annotation_width     = unit(c(0.8,0.4,0.6),"cm"),
  show_legend          = c(TRUE, TRUE, FALSE),
  gap = unit(2,"mm")
)

src_legend <- Legend(
  labels = names(src_cols),
  legend_gp = gpar(fill = src_cols),
  title = "Source",
  title_gp = gpar(fontsize=12, fontface="bold"),
  labels_gp = gpar(fontsize=11),
  grid_height = unit(4,"mm"),
  grid_width  = unit(4,"mm")
)

ht <- Heatmap(mat,
  name  = "Virulence",
  col   = c("0"="#F0F0F0","1"="#08519C"),
  left_annotation   = row_ann,
  show_row_names    = FALSE,
  show_row_dend     = FALSE,
  show_column_dend  = FALSE,
  column_names_gp   = gpar(fontsize=13, fontface="italic"),
  column_names_rot  = 40,
  column_names_side = "bottom",
  row_split         = ann$ST_split,
  row_gap           = unit(4,"mm"),
  row_title_gp      = gpar(fontsize=13, fontface="bold"),
  row_title_rot     = 0,
  column_title      = NULL,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  rect_gp = gpar(col="white", lwd=0.5),
  heatmap_legend_param = list(
    title="Present", labels=c("No","Yes"), at=c(0,1),
    legend_height=unit(2,"cm"),
    labels_gp=gpar(fontsize=12),
    title_gp=gpar(fontsize=12, fontface="bold")
  )
)

png(file.path(OUT,"Fig3_virulence_heatmap.png"),
    width=10, height=15, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,40,5),"mm"),
     annotation_legend_list=list(src_legend))
dev.off()

pdf(file.path(OUT,"Fig3_virulence_heatmap.pdf"), width=10, height=15)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,40,5),"mm"),
     annotation_legend_list=list(src_legend))
dev.off()

message("Fig3 saved")
