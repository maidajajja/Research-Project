suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

amr <- read.table("/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder_final229/amrfinder_all_final229.tsv",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE, quote="")
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

st_col <- grep("mlst__ST$", colnames(kleb), value=TRUE)[1]
kleb$ST <- kleb[[st_col]]

amr2 <- left_join(amr, kleb[,c("strain","ST")], by=c("Sample"="strain"))
amr2 <- amr2[!is.na(amr2$Type) & amr2$Type == "AMR", ]
amr2$ST[is.na(amr2$ST)] <- "Unknown"

meta_hh <- bind_rows(
  data.frame(key=as.character(meta[["Genome ID"]]), HH=meta[["Host_Health_Clean"]], stringsAsFactors=FALSE),
  data.frame(key=as.character(meta[["Assembly Accession"]]), HH=meta[["Host_Health_Clean"]], stringsAsFactors=FALSE),
  data.frame(key=paste0(as.character(meta[["SRA Accession"]]),"_assembled"), HH=meta[["Host_Health_Clean"]], stringsAsFactors=FALSE)
) %>% filter(!is.na(key) & key != "" & key != "NA_assembled" & key != "nan_assembled") %>%
  distinct(key, .keep_all=TRUE)

amr2$meta_key <- ifelse(grepl("^GC[AF]_", amr2$Sample),
  sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", amr2$Sample), amr2$Sample)
gca_rows <- meta_hh[grepl("^GC[AF]_", meta_hh$key), ]
gca_rows$key <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", gca_rows$key)
meta_hh <- bind_rows(meta_hh, gca_rows) %>% distinct(key, .keep_all=TRUE)
amr2 <- left_join(amr2, meta_hh, by=c("meta_key"="key"))
amr2$HH[is.na(amr2$HH)] <- "Unknown"

top_genes <- names(sort(table(amr2$Element.symbol), decreasing=TRUE))[1:25]
amr_filt <- amr2[amr2$Element.symbol %in% top_genes, ]

gene_class <- amr_filt %>%
  distinct(Element.symbol, Class) %>%
  mutate(Class_clean = case_when(
    Class == "AMINOGLYCOSIDE" ~ "Aminoglycoside",
    Class == "BETA-LACTAM" ~ "Beta-lactam",
    Class == "FOSFOMYCIN" ~ "Fosfomycin",
    Class == "MACROLIDE" ~ "Macrolide",
    grepl("QUINOLONE", Class) ~ "Quinolone",
    Class == "SULFONAMIDE" ~ "Sulfonamide",
    Class == "TETRACYCLINE" ~ "Tetracycline",
    TRUE ~ "Other"
  )) %>%
  group_by(Element.symbol) %>%
  slice(1) %>% ungroup()

pa <- amr_filt %>%
  distinct(Sample, Element.symbol) %>%
  mutate(present=1) %>%
  pivot_wider(names_from=Element.symbol, values_from=present, values_fill=0)

sample_meta <- amr_filt %>%
  distinct(Sample, ST, HH) %>%
  group_by(Sample) %>% slice(1) %>% ungroup()

pa2 <- left_join(pa, sample_meta, by="Sample")
pa2$ST[is.na(pa2$ST)] <- "Unknown"
pa2$HH[is.na(pa2$HH)] <- "Unknown"

# ST order by mean AMR gene count descending
top_sts <- c("ST258","ST11","ST512","ST65","ST23","ST86","ST29")
st_levels <- c(top_sts, "Other")

pa2$ST_group <- factor(
  ifelse(pa2$ST %in% top_sts, pa2$ST, "Other"),
  levels=st_levels)
pa2 <- pa2 %>% arrange(ST_group, desc(rowSums(across(all_of(top_genes)))))
split_vec <- factor(as.character(pa2$ST_group), levels=st_levels)

# Cap rows per ST to max 15 to prevent ST23/Other dominating
set.seed(42)
pa2_capped <- pa2 %>%
  group_by(ST_group) %>%
  slice_head(n=15) %>%
  ungroup()
split_capped <- factor(as.character(pa2_capped$ST_group), levels=st_levels)

cat("Capped ST counts:\n"); print(table(split_capped))

# Gene order by class then prevalence
gene_class_ordered <- gene_class %>%
  arrange(Class_clean, desc(sapply(Element.symbol,
    function(g) sum(pa2[[g]], na.rm=TRUE))))
gene_order <- gene_class_ordered$Element.symbol
gene_order <- gene_order[gene_order %in% colnames(pa2_capped)]

mat <- as.matrix(pa2_capped[, gene_order])
rownames(mat) <- pa2_capped$Sample

# Paul Tol muted palette for AMR classes - distinct from ST colours
class_colours <- c(
  "Aminoglycoside" = "#44AA99",
  "Beta-lactam"    = "#88CCEE",
  "Fosfomycin"     = "#117733",
  "Macrolide"      = "#999933",
  "Quinolone"      = "#CC6677",
  "Sulfonamide"    = "#AA4499",
  "Tetracycline"   = "#882255"
)

gene_classes <- gene_class_ordered$Class_clean[
  match(gene_order, gene_class_ordered$Element.symbol)]

col_ann <- HeatmapAnnotation(
  `AMR Class` = gene_classes,
  col = list(`AMR Class` = class_colours),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_name_side = "left",
  simple_anno_size = unit(7, "mm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

# ST colours - distinct from AMR class colours
# ST23 changed to orange-red to avoid clash with Macrolide yellow
st_cols <- setNames(c(
  "#332288",  # ST258 - dark purple
  "#117733",  # ST11  - dark green
  "#44AA99",  # ST512 - teal — wait this clashes with Aminoglycoside
  "#6699CC",  # ST65  - muted blue
  "#994455",  # ST23  - dark rose (NOT yellow - avoids Macrolide clash)
  "#CC6677",  # ST86  - muted pink-red
  "#AA4499",  # ST29  - muted purple-pink
  "#BBBBBB"   # Other - grey
), st_levels)

# Host Health - muted, distinct from both ST and AMR class colours
hh_cols <- c(
  "Liver abscess"            = "#D6604D",
  "Liver transplant"         = "#4393C3",
  "Other liver disease"      = "#74C476",
  "Non-liver or unspecified" = "#FDB863",
  "Unknown"                  = "#DDDDDD"
)
pa2_capped$HH <- factor(pa2_capped$HH, levels=names(hh_cols))

row_ann <- rowAnnotation(
  ST = split_capped,
  `Host Health` = pa2_capped$HH,
  col = list(
    ST = st_cols,
    `Host Health` = hh_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(7, "mm"),
  gap = unit(2, "mm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

ht <- Heatmap(
  mat,
  name = "AMR gene",
  col = c("0"="grey95", "1"="#5b5b5b"),
  rect_gp = gpar(col="white", lwd=1.0),
  top_annotation = col_ann,
  row_split = split_capped,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize=9, fontface="italic"),
  column_names_rot = 45,
  column_title = NULL,
  row_title_gp = gpar(fontsize=11, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(4, "mm"),
  left_annotation = row_ann,
  heatmap_legend_param = list(
    title = "AMR gene",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey95","#5b5b5b")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)),
  border = TRUE,
  border_gp = gpar(col="grey80", lwd=0.3),
  use_raster = FALSE)

png(file.path(OUT,"Fig2_AMR_heatmap_v4.png"),
    width=16, height=14, units="in", res=600)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,15),"mm"), background="white")
dev.off()

pdf(file.path(OUT,"Fig2_AMR_heatmap_v4.pdf"), width=16, height=14)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,15),"mm"))
dev.off()

message("Fig2 v4 saved")
