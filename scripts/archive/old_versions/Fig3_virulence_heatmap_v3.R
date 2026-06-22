suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
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
meta$Source <- dplyr::recode(meta$Source,
  "liver abscess"="Liver abscess","hepatic abscess"="Liver abscess",
  "liver abscess drainage fluid"="Drainage/Pus",
  "liver abscess puncture fluid"="Drainage/Pus",
  "Liver abscess puncture fluid"="Drainage/Pus",
  "abscess drainage"="Drainage/Pus","liver drainage"="Drainage/Pus",
  "liver"="Liver","blood"="Blood","Blood"="Blood")
meta$Health[is.na(meta$Health)|meta$Health==""] <- "Unknown"
meta <- meta[meta$Sample!="", c("Sample","Source","Health")]

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col  <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
kleb$ST  <- gsub("-.*","", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST,
                     Vir=kleb$Vir, stringsAsFactors=FALSE)

vir_cols <- c("klebsiella__ybst__Yersiniabactin",
              "klebsiella__cbst__Colibactin",
              "klebsiella__abst__Aerobactin",
              "klebsiella__smst__Salmochelin",
              "klebsiella__rmst__RmpADC",
              "klebsiella__rmpa2__rmpA2")
vir_cols <- vir_cols[vir_cols %in% colnames(kleb)]

vir_df <- kleb[, c("strain", vir_cols)]
mat <- matrix(0, nrow=nrow(vir_df), ncol=length(vir_cols),
              dimnames=list(vir_df$strain, vir_cols))
for(col in vir_cols){
  mat[,col] <- ifelse(is.na(vir_df[[col]]) | vir_df[[col]]=="-" |
                      vir_df[[col]]=="0" | vir_df[[col]]=="", 0, 1)
}
colnames(mat) <- c("Yersiniabactin","Colibactin","Aerobactin",
                   "Salmochelin","RmpADC","rmpA2")[1:ncol(mat)]

ann <- data.frame(Sample=rownames(mat), stringsAsFactors=FALSE) %>%
  left_join(st_map, by="Sample") %>%
  left_join(meta, by="Sample")
ann$ST[is.na(ann$ST)]         <- "Unknown"
ann$Source[is.na(ann$Source)] <- "Unknown"
ann$Health[is.na(ann$Health)] <- "Unknown"
ann$Vir[is.na(ann$Vir)]       <- 0

# ── ST ordering by mean virulence locus count ascending ───────────────────────
ann$gene_count <- rowSums(mat[match(ann$Sample, rownames(mat)), ])
st_mean <- ann %>%
  filter(ST != "Unknown") %>%
  group_by(ST) %>%
  summarise(mean_vir=mean(gene_count), .groups="drop") %>%
  arrange(mean_vir)

top_sts_by_count <- names(sort(table(ann$ST[ann$ST!="Unknown"]), decreasing=TRUE))[1:7]
st_mean_filtered <- st_mean %>% filter(ST %in% top_sts_by_count)
top_sts_asc <- st_mean_filtered$ST

st_pal  <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
st_cols <- setNames(c(st_pal,"grey80","grey90"),
                    c(top_sts_asc,"Other","Unknown"))
ann$ST_group <- ifelse(ann$ST %in% top_sts_asc, ann$ST,
                       ifelse(ann$ST=="Unknown","Unknown","Other"))
ann <- ann[!is.na(ann$ST_group) & ann$ST_group!="NA", ]
mat <- mat[rownames(mat) %in% ann$Sample, ]

ann$gene_count <- rowSums(mat[match(ann$Sample, rownames(mat)), ])
ann <- ann %>%
  arrange(factor(ST_group, levels=c(top_sts_asc,"Other","Unknown")), gene_count)
mat <- mat[match(ann$Sample, rownames(mat)), ]
ann$ST_split <- factor(ann$ST_group, levels=c(top_sts_asc,"Other","Unknown"))

src_levels <- c("Liver abscess","Blood","Drainage/Pus","Liver","Other","Unknown")
src_cols <- setNames(c("#D55E00","#CC79A7","#56B4E9","#009E73","#4D4D4D","grey85"),
                     src_levels)
ann$Src_clean <- ifelse(ann$Source %in% src_levels, ann$Source, "Other")

health_vals <- unique(ann$Health)
health_cols <- setNames(
  colorRampPalette(brewer.pal(8,"Accent"))(length(health_vals)),
  health_vals)

row_ann <- rowAnnotation(
  ST           = ann$ST_group,
  `Vir. Score` = ann$Vir,
  Source       = ann$Src_clean,
  col = list(
    ST           = st_cols,
    `Vir. Score` = colorRamp2(c(0,2,5), c("#FFF7BC","#FE9929","#CC4C02")),
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
  column_names_gp  = gpar(fontsize=10, fontface="bold"),
  column_names_rot = 45,
  row_split        = ann$ST_split,
  row_gap          = unit(3,"mm"),
  row_title_gp     = gpar(fontsize=9, fontface="bold"),
  row_title_rot    = 0,
  column_title     = NULL,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  rect_gp = gpar(col="white", lwd=0.5),
  heatmap_legend_param = list(
    title="Present", labels=c("No","Yes"),
    at=c(0,1), legend_height=unit(2,"cm")
  )
)

png(file.path(OUT,"Fig3_virulence_heatmap.png"),
    width=9, height=14, units="in", res=600, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,5),"mm"))
dev.off()
pdf(file.path(OUT,"Fig3_virulence_heatmap.pdf"), width=9, height=14)
draw(ht, merge_legend=TRUE)
dev.off()
message("Fig3 v3 saved at 600 DPI")
