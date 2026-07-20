suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# ── Load pan-genome matrix ────────────────────────────────────────────────────
pa <- read.table("/scratch/users/k22017808/KP_Research_Project/06_Pangenome_final229/gene_presence_absence.Rtab",
                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)
mat <- as.matrix(pa)
n   <- ncol(mat)
freq <- rowSums(mat) / n

cat("Total genes:", nrow(mat), "Isolates:", n, "\n")

# Gene categories
gene_cat <- ifelse(freq >= 0.99, "Core",
            ifelse(freq >= 0.95, "Soft core",
            ifelse(freq >= 0.15, "Shell", "Cloud")))

# ── Load metadata ─────────────────────────────────────────────────────────────
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col  <- grep("mlst__ST$", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST  <- sub("^ST","", gsub("-.*","", kleb[[st_col]]))
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))

# Match sample metadata
top_sts <- c("23","11","86","258","65","29","512")

# Build fast lookup from genomes.csv
# Key: short accession (GCA_xxx.1) -> metadata
meta$acc_short <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1",
                       trimws(meta[["Assembly Accession"]]))
meta$sra_key <- paste0(trimws(meta[["SRA Accession"]]),"_assembled")
meta$gid_key <- as.character(meta[["Genome ID"]])

get_meta <- function(sample) {
  # Try Kleborate strain match first
  k <- kleb[kleb$strain == sample, ]
  if(nrow(k) > 0) {
    st  <- k$ST[1]
    vir <- k$Vir[1]
    res <- k$Res[1]
  } else {
    st <- "Unknown"; vir <- NA; res <- NA
  }
  
  # Extract short accession from full sample name
  # e.g. GCA_019165765.1_ASM1916576v1_genomic -> GCA_019165765.1
  sample_short <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", sample)
  
  # Match via short accession, SRA, or Genome ID
  m <- meta[meta$acc_short == sample_short |
            meta$sra_key == sample |
            meta$gid_key == sample, ]
  
  hh   <- if(nrow(m)>0) m$Host_Health_Clean[1] else "Unknown"
  cont <- if(nrow(m)>0) m$Continent[1] else "Unknown"
  if(is.na(hh)||hh==""|hh=="nan") hh <- "Unknown"
  if(is.na(cont)||cont==""|cont=="nan") cont <- "Unknown"
  
  list(ST=st, Vir=vir, Res=res, HH=hh, Cont=cont)
}

samples <- colnames(mat)
meta_list <- lapply(samples, get_meta)
sample_st   <- sapply(meta_list, `[[`, "ST")
sample_vir  <- sapply(meta_list, `[[`, "Vir")
sample_res  <- sapply(meta_list, `[[`, "Res")
sample_hh   <- sapply(meta_list, `[[`, "HH")
sample_cont <- sapply(meta_list, `[[`, "Cont")

sample_st_group <- ifelse(sample_st %in% top_sts,
                          paste0("ST", sample_st), "Other")

cat("ST distribution:\n"); print(table(sample_st_group))
cat("HH distribution:\n"); print(table(sample_hh))

# Order samples by ST then virulence score
st_order <- c(paste0("ST", top_sts), "Other")
st_order_plot <- st_order  # No clustering so order is preserved top-to-bottom
sample_order <- order(
  factor(sample_st_group, levels=st_order),
  -as.numeric(sample_vir))
# ST23 will be first as top_sts[1]='23'
mat_ord      <- mat[, sample_order]
st_ord       <- sample_st_group[sample_order]
vir_ord      <- as.numeric(sample_vir)[sample_order]
res_ord      <- as.numeric(sample_res)[sample_order]
hh_ord       <- sample_hh[sample_order]
cont_ord     <- sample_cont[sample_order]

# Paul Tol muted ST colours
st_cols <- setNames(c(
  "#774411","#225522","#CC99BB","#332288",
  "#EE8866","#AAAA00","#77AADD","#BBBBBB"),
  st_order)

hh_cols <- c(
  "Liver abscess"      = "#994455",
  "Liver transplant"   = "#4477AA",
  "Other liver disease"= "#228833",
  "Unknown"            = "#DDDDDD")

cont_cols <- c(
  "Asia"         = "#DDAA33",
  "Europe"       = "#BB5566",
  "North America"= "#004488",
  "Unknown"      = "#BBBBBB")

# Row annotation (isolates)
col_ann <- rowAnnotation(
  ST = st_ord,
  `Host Health` = hh_ord,
  Continent = cont_ord,
  col = list(
    ST = st_cols,
    `Host Health` = hh_cols,
    Continent = cont_cols),
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize=16, fontface="bold"),
  gap = unit(1,"mm"),
  annotation_legend_param = list(
    ST = list(title_gp=gpar(fontsize=16,fontface="bold"),
              labels_gp=gpar(fontsize=16)),
    `Host Health` = list(title_gp=gpar(fontsize=16,fontface="bold"),
                         labels_gp=gpar(fontsize=16)),
    Continent = list(title_gp=gpar(fontsize=16,fontface="bold"),
                     labels_gp=gpar(fontsize=16))))

# ── Panel A: Core + Soft core ────────────────────────────────────────────────
core_sc_idx <- which(gene_cat %in% c("Core","Soft core"))
mat_A <- mat_ord[core_sc_idx, ]
cat_A <- gene_cat[core_sc_idx]

# Order genes by frequency descending
gene_order_A <- order(-rowSums(mat_A))
mat_A <- mat_A[gene_order_A, ]
cat_A <- cat_A[gene_order_A]

# Top annotation — gene category + prevalence
top_ann_A <- HeatmapAnnotation(
  Category = cat_A,
  Prevalence = anno_barplot(rowSums(mat_A)/ncol(mat_A),
                            height=unit(3.5,"cm"),
                            ylim=c(-0.05,1), gp=gpar(fill="grey50", col=NA), axis_param=list(side="left", at=c(0,0.5,1), labels=c("0","0.5","1"), gp=gpar(fontsize=16))),
  col = list(Category = c("Core"="#1A5276","Soft core"="#2E86C1","Shell"="#E67E22")),
  annotation_legend_param = list(Category = list(at=c("Core","Soft core","Shell"), labels=c("Core","Soft core","Shell"), title_gp=gpar(fontsize=16,fontface="bold"), labels_gp=gpar(fontsize=16))),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize=12, fontface="bold"),
  show_annotation_name = c(Category = TRUE, Prevalence = TRUE),
  gap = unit(1,"mm"))

htA <- Heatmap(t(mat_A),
  name = "Gene presence",
  col = c("0"="grey97","1"="#2C3E50"),
  top_annotation = top_ann_A,
  left_annotation = col_ann,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  column_split = cat_A,
  column_gap = unit(2,"mm"),
  column_title_gp = gpar(fontsize=16, fontface="bold"),
  rect_gp = gpar(col=NA),
  heatmap_legend_param = list(
    title="Gene presence",
    labels=c("Absent","Present"),
    at=c(0,1),
    legend_gp=gpar(fill=c("grey97","#2C3E50")),
    title_gp=gpar(fontsize=16,fontface="bold"),
    labels_gp=gpar(fontsize=16)),
  row_split = factor(st_ord, levels=st_order_plot),
  row_title = NULL,
  row_gap = unit(2,"mm"),
  row_title_gp = gpar(fontsize=16, fontface="bold"),
  row_title_rot = 0,
  border = TRUE,
  border_gp = gpar(col="grey85", lwd=0.3),
  column_title = NULL,
  column_title_side = "bottom",
  use_raster = TRUE,
  raster_quality = 2)

# ── Panel B: Shell genes ──────────────────────────────────────────────────────
shell_idx <- which(gene_cat == "Shell")
mat_B <- mat_ord[shell_idx, ]
cat_B <- gene_cat[shell_idx]

gene_order_B <- order(-rowSums(mat_B))
mat_B <- mat_B[gene_order_B, ]
cat_B <- cat_B[gene_order_B]

top_ann_B <- HeatmapAnnotation(
  Category = cat_B,
  `Prevalence ` = anno_barplot(rowSums(mat_B)/ncol(mat_B),
                            height=unit(3.5,"cm"),
                            ylim=c(-0.05,1), gp=gpar(fill="grey50", col=NA), axis_param=list(side="left", at=c(0,0.5,1), labels=c("","",""), gp=gpar(fontsize=16))),
  col = list(Category = c("Shell"="#E67E22")),
  annotation_name_side = "right",
  show_annotation_name = c(Category = FALSE, `Prevalence ` = FALSE),
  annotation_name_gp = gpar(fontsize=16, fontface="bold"),
  gap = unit(1,"mm"),
  show_legend = FALSE)

htB <- Heatmap(t(mat_B),
  name = "Gene presence (shell)",
  col = c("0"="grey97","1"="#2C3E50"),
  top_annotation = top_ann_B,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  rect_gp = gpar(col=NA),
  show_heatmap_legend = FALSE,
  row_split = factor(st_ord, levels=st_order_plot),
  row_title = NULL,
  row_gap = unit(2,"mm"),
  row_title_gp = gpar(fontsize=16, fontface="bold"),
  row_title_rot = 0,
  border = TRUE,
  border_gp = gpar(col="grey85", lwd=0.3),
  column_title = NULL,
  column_title_side = "bottom",
  use_raster = TRUE,
  raster_quality = 2)

# ── Save ──────────────────────────────────────────────────────────────────────
cat("Drawing heatmap...\n")
png(file.path(OUT,"Fig6c_pangenome_heatmap_v4.png"),
    width=20, height=12, units="in", res=300, bg="white")
draw(htA + htB,
     merge_legend=TRUE,
     padding=unit(c(25,25,25,60),"mm"),
     column_title=NULL,
     column_title_gp=gpar(fontsize=16, fontface="bold"))
dev.off()

pdf(file.path(OUT,"Fig6c_pangenome_heatmap_v4.pdf"),
    width=20, height=12)
draw(htA + htB,
     merge_legend=TRUE,
     padding=unit(c(25,25,25,60),"mm"),
     column_title=NULL,
     column_title_gp=gpar(fontsize=16, fontface="bold"))
dev.off()

message("Fig6c v3 saved — two-panel Ellis-style heatmap")
