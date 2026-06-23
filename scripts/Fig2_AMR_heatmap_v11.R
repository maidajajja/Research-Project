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

# Point mutations from AMRFinder method
point_mutation_genes <- amr2 %>%
  filter(Method %in% c("POINTX","POINTN")) %>%
  distinct(Element.symbol) %>% pull()

# Intrinsic chromosomal genes in Kp - manually curated
# ompK36_D135DGD is a POINT MUTATION in a porin - belongs in Quinolone/Beta-lactam not intrinsic
intrinsic_kp <- c("fosA","oqxA","oqxB","oqxA11","oqxB19","blaSHV-1","blaSHV-11")

top_genes <- names(sort(table(amr2$Element.symbol), decreasing=TRUE))[1:25]
amr_filt <- amr2[amr2$Element.symbol %in% top_genes, ]

# Gene class mapping - ompK36 now correctly in Beta-lactam (porin affecting carbapenem entry)
gene_class <- amr_filt %>%
  distinct(Element.symbol, Class) %>%
  mutate(Class_clean = case_when(
    Element.symbol %in% intrinsic_kp          ~ "Intrinsic*",
    Element.symbol == "ompK36_D135DGD"        ~ "Beta-lactam\u2020",  # porin mutation
    Class == "AMINOGLYCOSIDE"                  ~ "Aminoglycoside",
    Class == "BETA-LACTAM"                     ~ "Beta-lactam",
    Class == "MACROLIDE"                       ~ "Macrolide",
    grepl("QUINOLONE", Class)                  ~ "Quinolone",
    Class == "SULFONAMIDE"                     ~ "Sulfonamide",
    Class == "TETRACYCLINE"                    ~ "Tetracycline",
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

top_sts <- c("ST258","ST11","ST512","ST65","ST23","ST86","ST29")
pa2 <- pa2 %>% filter(ST %in% top_sts)
pa2$ST_group <- factor(pa2$ST, levels=top_sts)

acquired_in_mat <- setdiff(top_genes, intrinsic_kp)
acquired_in_mat <- acquired_in_mat[acquired_in_mat %in% colnames(pa2)]
pa2 <- pa2 %>% arrange(ST_group, desc(rowSums(across(all_of(acquired_in_mat)))))
split_vec <- factor(as.character(pa2$ST_group), levels=top_sts)

# ST labels with n
st_counts <- table(split_vec)
st_labels <- paste0(top_sts, " (n=", st_counts[top_sts], ")")
names(st_labels) <- top_sts

# Gene ordering
class_levels <- c("Aminoglycoside","Beta-lactam","Beta-lactam\u2020",
                  "Macrolide","Quinolone","Sulfonamide","Tetracycline","Intrinsic*")

gene_class_ordered <- gene_class %>%
  mutate(Class_clean = factor(Class_clean, levels=class_levels)) %>%
  arrange(Class_clean, desc(sapply(Element.symbol,
    function(g) if(g %in% colnames(pa2)) sum(pa2[[g]], na.rm=TRUE) else 0)))
gene_order <- gene_class_ordered$Element.symbol
gene_order <- gene_order[gene_order %in% colnames(pa2)]

# Column split - merge Beta-lactam and Beta-lactam† for display
col_split_raw <- gene_class_ordered$Class_clean[match(gene_order, gene_class_ordered$Element.symbol)]
col_split <- factor(as.character(col_split_raw), levels=class_levels)

mat <- as.matrix(pa2[, gene_order])
rownames(mat) <- pa2$Sample

# Gene labels with markers
gene_labels <- ifelse(gene_order %in% point_mutation_genes,
                      paste0(gene_order, "\u2020"),
               ifelse(gene_order %in% intrinsic_kp,
                      paste0(gene_order, "*"),
                      gene_order))

# Prevalence across all 129 isolates shown
gene_prev <- colSums(mat) / nrow(mat) * 100

# Acquired genes per isolate
acquired_cols <- gene_order[!(gene_order %in% intrinsic_kp)]
total_acquired <- rowSums(mat[, acquired_cols, drop=FALSE])

# Colours - no Fosfomycin class needed anymore
class_colours <- c(
  "Aminoglycoside"   = "#44AA99",
  "Beta-lactam"      = "#6699CC",
  "Beta-lactam\u2020"= "#336699",  # darker blue for porin mutation
  "Macrolide"        = "#EEBB88",
  "Quinolone"        = "#CC6677",
  "Sulfonamide"      = "#AA4499",
  "Tetracycline"     = "#882255",
  "Intrinsic*"       = "#BBBBBB"
)

# Prevalence bar colours matching class
prev_cols <- sapply(as.character(col_split), function(x) class_colours[x])

# Top annotation - prevalence bar only, no overlapping label
top_ann <- HeatmapAnnotation(
  `Prevalence\n(% isolates)` = anno_barplot(
    gene_prev,
    gp = gpar(fill=prev_cols, col=NA),
    height = unit(20, "mm"),
    ylim = c(0, 100),
    axis_param = list(
      side = "left",
      gp = gpar(fontsize=7),
      at = c(0,25,50,75,100),
      labels = c("0","25","50","75","100%"))),
  `AMR Class` = col_split,
  col = list(`AMR Class` = class_colours),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  simple_anno_size = unit(5, "mm"),
  gap = unit(2, "mm"),
  annotation_legend_param = list(
    `AMR Class` = list(
      title = "AMR Class",
      title_gp = gpar(fontsize=9, fontface="bold"),
      labels_gp = gpar(fontsize=8))))

# Right annotation - acquired genes bar
right_ann <- rowAnnotation(
  `Acquired genes` = anno_barplot(
    total_acquired,
    gp = gpar(fill="#5b5b5b", col=NA),
    width = unit(25, "mm"),
    axis_param = list(
      side = "top",
      gp = gpar(fontsize=7))),
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  annotation_name_side = "top",
  annotation_name_rot = 0)

# Left annotation
st_cols <- setNames(c(
  "#332288","#225522","#77AADD","#EE8866",
  "#774411","#CC99BB","#AAAA00"), top_sts)

hh_cols <- c(
  "Liver abscess"            = "#994455",
  "Liver transplant"         = "#4477AA",
  "Other liver disease"      = "#228833",
  "Non-liver or unspecified" = "#DDAA33",
  "Unknown"                  = "#DDDDDD")
pa2$HH <- factor(pa2$HH, levels=names(hh_cols))

row_ann <- rowAnnotation(
  ST = split_vec,
  `Host Health` = pa2$HH,
  col = list(ST=st_cols, `Host Health`=hh_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(6, "mm"),
  gap = unit(2, "mm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

n_rows <- nrow(mat)
fig_height <- max(14, n_rows * 0.09 + 6)

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
  column_labels = gene_labels,
  column_names_side = "top",
  column_names_gp = gpar(fontsize=8, fontface="italic"),
  column_names_rot = 45,
  column_title = NULL,
  row_title = st_labels[top_sts],
  row_title_gp = gpar(fontsize=10, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(4, "mm"),
  height = unit(n_rows * 2.5, "mm"),
  left_annotation = row_ann,
  right_annotation = right_ann,
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

png(file.path(OUT,"Fig2_AMR_heatmap_v11.png"),
    width=18, height=fig_height, units="in", res=300)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,15),"mm"), background="white")
dev.off()

pdf(file.path(OUT,"Fig2_AMR_heatmap_v11.pdf"), width=18, height=fig_height)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,15),"mm"))
dev.off()

message(paste0("Fig2 v11 saved - ", n_rows, " isolates"))
