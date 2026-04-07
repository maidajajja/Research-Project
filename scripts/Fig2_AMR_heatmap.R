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
meta$Country <- meta[["Isolation Country"]]
meta$Source  <- meta[["Isolation Source"]]
meta$Country[is.na(meta$Country)|meta$Country==""] <- "Unknown"
meta$Source[is.na(meta$Source)|meta$Source==""]   <- "Unknown"
meta <- meta[meta$Sample!="", c("Sample","Country","Source")]

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col  <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST  <- gsub("-.*","", kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))
st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST,
                     Vir=kleb$Vir, Res=kleb$Res, stringsAsFactors=FALSE)

amr_dir <- "/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder"
files <- list.files(amr_dir, pattern="_amrfinder\\.tsv$", full.names=TRUE)
amr_list <- lapply(files, function(f) {
  sample <- gsub("_amrfinder\\.tsv$","", basename(f))
  df <- tryCatch(read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                             quote="", check.names=FALSE), error=function(e) NULL)
  if (is.null(df)||nrow(df)==0||!"Element symbol" %in% colnames(df)) return(NULL)
  data.frame(Sample=sample, Gene=df[["Element symbol"]], Class=df[["Class"]])
})
amr_all <- do.call(rbind, Filter(Negate(is.null), amr_list))
amr_all$Class <- dplyr::recode(amr_all$Class,
  "NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE" = "MULTI-CLASS",
  "PHENICOL/QUINOLONE" = "PHENICOL/QUIN")

gene_prev <- amr_all %>% group_by(Gene) %>%
  summarise(n=n_distinct(Sample)) %>% arrange(desc(n)) %>% slice_head(n=25)
amr_filt <- amr_all %>% filter(Gene %in% gene_prev$Gene)

amr_wide <- amr_filt %>% mutate(present=1) %>%
  pivot_wider(id_cols=Sample, names_from=Gene, values_from=present,
              values_fill=0, values_fn=max) %>% as.data.frame()
rownames(amr_wide) <- amr_wide$Sample
mat <- as.matrix(amr_wide[,-1])

ann <- data.frame(Sample=rownames(mat), stringsAsFactors=FALSE) %>%
  left_join(st_map, by="Sample") %>%
  left_join(meta, by="Sample")
ann$ST[is.na(ann$ST)]           <- "Unknown"
ann$Country[is.na(ann$Country)] <- "Unknown"
ann$Source[is.na(ann$Source)]   <- "Unknown"
ann$Vir[is.na(ann$Vir)]         <- 0
ann$Res[is.na(ann$Res)]         <- 0

# ST colours — wide bar like supervisor
top_sts <- names(sort(table(ann$ST[ann$ST!="Unknown"]), decreasing=TRUE))[1:7]
st_pal  <- brewer.pal(7,"Set1")
st_cols <- setNames(c(st_pal,"grey80","grey90"), c(top_sts,"Other","Unknown"))
ann$ST_group <- ifelse(ann$ST %in% top_sts, ann$ST,
                       ifelse(ann$ST=="Unknown","Unknown","Other"))

# Row split by top STs
ann$ST_split <- factor(ann$ST_group,
                       levels=c(top_sts[1:5],"Other","Unknown"))

# Left annotation — ST as wide bar, virulence and resistance as continuous
row_ann <- rowAnnotation(
  ST            = ann$ST_group,
  `Vir. Score`  = ann$Vir,
  `Res. Score`  = ann$Res,
  col = list(
    ST           = st_cols,
    `Vir. Score` = colorRamp2(c(0,2,5), c("#F7FBFF","#6BAED6","#08306B")),
    `Res. Score` = colorRamp2(c(0,1,3), c("#FFF5F0","#FB6A4A","#67000D"))
  ),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_width   = unit(c(0.8,0.4,0.4),"cm"),
  gap = unit(1,"mm")
)

ht <- Heatmap(mat,
  name  = "AMR gene",
  col   = c("0"="#F0F0F0","1"="#1a1a1a"),
  left_annotation  = row_ann,
  show_row_names   = FALSE,
  show_column_dend = TRUE,
  column_names_gp  = gpar(fontsize=7.5, fontface="italic"),
  column_names_rot = 45,
  row_split        = ann$ST_split,
  row_gap          = unit(3,"mm"),
  row_title_gp     = gpar(fontsize=9, fontface="bold"),
  row_title_rot    = 0,
  column_title     = expression(italic("K. pneumoniae")~"AMR Gene Presence/Absence (n = 234)"),
  column_title_gp  = gpar(fontsize=12, fontface="bold"),
  clustering_distance_rows    = "binary",
  clustering_distance_columns = "binary",
  clustering_method_rows      = "ward.D2",
  clustering_method_columns   = "ward.D2",
  rect_gp = gpar(col="white", lwd=0.3),
  heatmap_legend_param = list(
    title="Present",
    labels=c("No","Yes"),
    at=c(0,1),
    legend_height=unit(2,"cm")
  )
)

png(file.path(OUT,"Fig2_AMR_heatmap.png"),
    width=13, height=14, units="in", res=300, bg="white")
draw(ht, merge_legend=TRUE,
     padding=unit(c(5,5,5,5),"mm"))
dev.off()

pdf(file.path(OUT,"Fig2_AMR_heatmap.pdf"),
    width=13, height=14)
draw(ht, merge_legend=TRUE)
dev.off()

message("Fig2 saved")
