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

# Manually define intrinsic/chromosomal genes in Kp based on published literature
# fosA, oqxA, oqxB, oqxA11, oqxB19, blaSHV-1, blaSHV-11 are chromosomally encoded in Kp
# gyrA/parC point mutations are acquired mutations not genes
# ompK36 is a porin modification
intrinsic_kp <- c("fosA", "oqxA", "oqxB", "oqxA11", "oqxB19",
                   "blaSHV-1", "blaSHV-11", "ompK36_D135DGD")

top_genes <- names(sort(table(amr2$Element.symbol), decreasing=TRUE))[1:25]
amr_filt <- amr2[amr2$Element.symbol %in% top_genes, ]

acquired_genes <- setdiff(top_genes, intrinsic_kp)
cat("Intrinsic genes:", intersect(top_genes, intrinsic_kp), "\n")
cat("Acquired genes:", acquired_genes, "\n")

gene_class <- amr_filt %>%
  distinct(Element.symbol, Class) %>%
  mutate(Class_clean = case_when(
    Class == "AMINOGLYCOSIDE" ~ "Aminoglycoside",
    Class == "BETA-LACTAM"    ~ "Beta-lactam",
    Class == "FOSFOMYCIN"     ~ "Fosfomycin*",
    Class == "MACROLIDE"      ~ "Macrolide",
    grepl("QUINOLONE", Class) ~ "Quinolone",
    Class == "SULFONAMIDE"    ~ "Sulfonamide",
    Class == "TETRACYCLINE"   ~ "Tetracycline",
    TRUE ~ "Other"
  )) %>%
  mutate(Class_clean = ifelse(Element.symbol %in% intrinsic_kp,
                               "Intrinsic*", Class_clean)) %>%
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

# 7 key STs only
top_sts <- c("ST258","ST11","ST512","ST65","ST23","ST86","ST29")
pa2 <- pa2 %>% filter(ST %in% top_sts)
pa2$ST_group <- factor(pa2$ST, levels=top_sts)

# Order by acquired gene count (not total - so intrinsic don't inflate order)
acquired_in_mat <- acquired_genes[acquired_genes %in% colnames(pa2)]
pa2 <- pa2 %>% arrange(ST_group, desc(rowSums(across(all_of(acquired_in_mat)))))
split_vec <- factor(as.character(pa2$ST_group), levels=top_sts)

cat("ST counts:\n"); print(table(split_vec))

# Gene order: acquired genes by class, then intrinsic at right
class_levels <- c("Aminoglycoside","Beta-lactam","Macrolide","Quinolone",
                  "Sulfonamide","Tetracycline","Fosfomycin*","Intrinsic*")

gene_class_ordered <- gene_class %>%
  mutate(Class_clean = factor(Class_clean, levels=class_levels)) %>%
  arrange(Class_clean, desc(sapply(Element.symbol,
    function(g) if(g %in% colnames(pa2)) sum(pa2[[g]], na.rm=TRUE) else 0)))
gene_order <- gene_class_ordered$Element.symbol
gene_order <- gene_order[gene_order %in% colnames(pa2)]

col_split <- factor(
  gene_class_ordered$Class_clean[match(gene_order, gene_class_ordered$Element.symbol)],
  levels=class_levels)

mat <- as.matrix(pa2[, gene_order])
rownames(mat) <- pa2$Sample

cat("Col split:\n"); print(table(col_split))

# Paul Tol muted palette - all colourblind friendly and distinct
class_colours <- c(
  "Aminoglycoside" = "#44AA99",  # teal
  "Beta-lactam"    = "#6699CC",  # medium blue - distinct from teal
  "Macrolide"      = "#EEBB88",  # muted orange
  "Quinolone"      = "#CC6677",  # muted rose
  "Sulfonamide"    = "#AA4499",  # muted purple
  "Tetracycline"   = "#882255",  # dark magenta
  "Fosfomycin*"    = "#117733",  # dark green
  "Intrinsic*"     = "#BBBBBB"   # grey - visually separated
)

top_ann <- HeatmapAnnotation(
  `AMR Class` = col_split,
  col = list(`AMR Class` = class_colours),
  which = "column",
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  simple_anno_size = unit(5, "mm"),
  annotation_legend_param = list(
    title = "AMR Class",
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

# ST colours - all distinct and visible
st_cols <- setNames(c(
  "#332288",  # ST258 - dark indigo
  "#225522",  # ST11  - dark forest green
  "#77AADD",  # ST512 - medium blue
  "#EE8866",  # ST65  - muted orange
  "#774411",  # ST23  - dark brown
  "#CC99BB",  # ST86  - muted mauve
  "#AAAA00"   # ST29  - olive (visible, distinct)
), top_sts)

hh_cols <- c(
  "Liver abscess"            = "#994455",
  "Liver transplant"         = "#4477AA",
  "Other liver disease"      = "#228833",
  "Non-liver or unspecified" = "#DDAA33",
  "Unknown"                  = "#DDDDDD"
)
pa2$HH <- factor(pa2$HH, levels=names(hh_cols))

row_ann <- rowAnnotation(
  ST = split_vec,
  `Host Health` = pa2$HH,
  col = list(ST = st_cols, `Host Health` = hh_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(6, "mm"),
  gap = unit(2, "mm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

n_rows <- nrow(mat)
fig_height <- max(12, n_rows * 0.09 + 4)

ht <- Heatmap(
  mat,
  name = "AMR gene",
  col = c("0"="grey96", "1"="#5b5b5b"),
  rect_gp = gpar(col="white", lwd=0.8),
  top_annotation = top_ann,
  row_split = split_vec,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  column_split = col_split,
  column_gap = unit(3, "mm"),
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  show_row_names = FALSE,
  column_names_side = "top",
  column_names_gp = gpar(fontsize=8, fontface="italic"),
  column_names_rot = 45,
  column_title = NULL,
  row_title_gp = gpar(fontsize=11, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(4, "mm"),
  height = unit(n_rows * 2.5, "mm"),
  left_annotation = row_ann,
  heatmap_legend_param = list(
    title = "AMR gene",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey96","#5b5b5b")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)),
  border = TRUE,
  border_gp = gpar(col="grey85", lwd=0.3),
  use_raster = TRUE)

png(file.path(OUT,"Fig2_AMR_heatmap_v9.png"),
    width=16, height=fig_height, units="in", res=300)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,15),"mm"), background="white")
dev.off()

pdf(file.path(OUT,"Fig2_AMR_heatmap_v9.pdf"), width=16, height=fig_height)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,15),"mm"))
dev.off()

message(paste0("Fig2 v9 saved - ", n_rows, " isolates, intrinsic separated manually"))
