suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

# Extract columns
st_col   <- grep("mlst__ST$", colnames(kleb), value=TRUE)[1]
vir_col  <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
ybt_col  <- "klebsiella__ybst__Yersiniabactin"
clb_col  <- grep("cbst__Colibactin$", colnames(kleb), value=TRUE)[1]
iuc_col  <- grep("abst__Aerobactin$", colnames(kleb), value=TRUE)[1]
iro_col  <- grep("smst__Salmochelin$", colnames(kleb), value=TRUE)[1]
rmp_col  <- grep("rmst__RmpADC$", colnames(kleb), value=TRUE)[1]
rmp2_col <- grep("rmpa2__rmpA2$", colnames(kleb), value=TRUE)[1]

kleb$ST        <- kleb[[st_col]]
kleb$vir_score <- suppressWarnings(as.numeric(kleb[[vir_col]]))

# Binary presence/absence
kleb$Yersiniabactin <- ifelse(kleb[[ybt_col]] == "-", 0, 1)
kleb$Aerobactin     <- ifelse(kleb[[iuc_col]] == "-", 0, 1)
kleb$Salmochelin    <- ifelse(kleb[[iro_col]] == "-", 0, 1)
kleb$RmpADC         <- ifelse(kleb[[rmp_col]] == "-", 0, 1)
kleb$rmpA2          <- ifelse(kleb[[rmp2_col]] == "-", 0, 1)
kleb$Colibactin     <- ifelse(kleb[[clb_col]] == "-", 0, 1)

# Variant types - biologically meaningful
kleb$iuc_type <- case_when(
  kleb[[iuc_col]] == "-"           ~ "Absent",
  grepl("iuc 1", kleb[[iuc_col]])  ~ "iuc1",
  grepl("iuc 2", kleb[[iuc_col]])  ~ "iuc2",
  grepl("iuc 3", kleb[[iuc_col]])  ~ "iuc3",
  TRUE ~ "Other")

kleb$rmp_type <- case_when(
  kleb[[rmp_col]] == "-"           ~ "Absent",
  grepl("rmp 1", kleb[[rmp_col]])  ~ "rmp1 (KpVP-1)",
  grepl("rmp 2", kleb[[rmp_col]])  ~ "rmp2 (KpVP-2)",
  grepl("rmp 3", kleb[[rmp_col]])  ~ "rmp3 (ICEKp1)",
  TRUE ~ "Other")

kleb$ybt_type <- case_when(
  kleb[[ybt_col]] == "-"                ~ "Absent",
  grepl("plasmid", kleb[[ybt_col]])     ~ "Plasmid-encoded",
  grepl("ICEKp", kleb[[ybt_col]])       ~ "ICEKp-encoded",
  TRUE ~ "Other")

kleb$rmpA2_present <- ifelse(kleb[[rmp2_col]] == "-", "Absent",
                      ifelse(grepl("\\*", kleb[[rmp2_col]]), "Truncated",
                             "Intact"))

# Host Health lookup
meta_hh <- bind_rows(
  data.frame(key=as.character(meta[["Genome ID"]]), HH=meta[["Host_Health_Clean"]], stringsAsFactors=FALSE),
  data.frame(key=as.character(meta[["Assembly Accession"]]), HH=meta[["Host_Health_Clean"]], stringsAsFactors=FALSE),
  data.frame(key=paste0(as.character(meta[["SRA Accession"]]),"_assembled"), HH=meta[["Host_Health_Clean"]], stringsAsFactors=FALSE)
) %>% filter(!is.na(key) & key != "" & key != "NA_assembled" & key != "nan_assembled") %>%
  distinct(key, .keep_all=TRUE)

kleb$meta_key <- ifelse(grepl("^GC[AF]_", kleb$strain),
  sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", kleb$strain), kleb$strain)
gca_rows <- meta_hh[grepl("^GC[AF]_", meta_hh$key), ]
gca_rows$key <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", gca_rows$key)
meta_hh <- bind_rows(meta_hh, gca_rows) %>% distinct(key, .keep_all=TRUE)
kleb <- left_join(kleb, meta_hh, by=c("meta_key"="key"))
kleb$HH[is.na(kleb$HH)] <- "Unknown"

# 7 key STs ordered by DECREASING virulence
top_sts <- c("ST23","ST11","ST65","ST86","ST29","ST258","ST512")
kd <- kleb %>% filter(ST %in% top_sts)
kd$ST_group <- factor(kd$ST, levels=top_sts)
kd <- kd %>% arrange(ST_group, desc(vir_score))
split_vec <- factor(as.character(kd$ST_group), levels=top_sts)

st_counts <- table(split_vec)
st_labels <- paste0(top_sts, " (n=", st_counts[top_sts], ")")
names(st_labels) <- top_sts

# Matrix - 6 loci columns
loci <- c("Yersiniabactin","Aerobactin","Salmochelin","RmpADC","rmpA2","Colibactin")
mat <- as.matrix(kd[, loci])
rownames(mat) <- kd$strain

cat("Loci prevalence (%):\n")
print(round(colSums(mat)/nrow(mat)*100, 1))
cat("ST distribution:\n")
print(table(split_vec))

# Prevalence per locus
loci_prev <- colSums(mat) / nrow(mat) * 100

# Virulence scores for right bar
vir_scores <- kd$vir_score

# Paul Tol muted palette - one colour per locus
# Each locus gets its own distinct colour for the prevalence bar
loci_colours <- c(
  "Yersiniabactin" = "#44AA99",  # teal
  "Aerobactin"     = "#CC6677",  # muted rose
  "Salmochelin"    = "#882255",  # dark magenta
  "RmpADC"         = "#6699CC",  # medium blue
  "rmpA2"          = "#332288",  # dark indigo
  "Colibactin"     = "#997700"   # dark gold
)

# Cell colours per locus - each column gets its own colour when present
# Create a multi-layer heatmap approach
# Use a single matrix with per-cell colours

# Build colour matrix
col_mat <- matrix("grey96", nrow=nrow(mat), ncol=ncol(mat))
colnames(col_mat) <- loci
for(l in loci) {
  col_mat[mat[,l]==1, l] <- loci_colours[l]
}

# Top annotation - prevalence bars coloured by locus
top_ann <- HeatmapAnnotation(
  `Prevalence\n(%)` = anno_barplot(
    loci_prev,
    gp = gpar(fill=loci_colours[loci], col=NA),
    height = unit(25, "mm"),
    ylim = c(0,100),
    axis_param = list(
      side = "left",
      gp = gpar(fontsize=7),
      at = c(0,50,100),
      labels = c("0","50","100"))),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  show_legend = FALSE)

# Right annotation - virulence score as gradient bar
vir_col_fun <- colorRamp2(c(0,1,2,3,4,5),
  c("grey95","#FFEE99","#EEBB88","#CC6677","#994455","#661122"))

right_ann <- rowAnnotation(
  `Virulence\nscore` = anno_barplot(
    vir_scores,
    gp = gpar(fill=vir_col_fun(vir_scores), col=NA),
    width = unit(25, "mm"),
    ylim = c(0,5),
    axis_param = list(
      side = "bottom",
      gp = gpar(fontsize=7),
      at = c(0,1,2,3,4,5))),
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  width = unit(35, "mm"))

# ST and Host Health colours
st_cols <- setNames(c(
  "#774411",  # ST23
  "#225522",  # ST11
  "#EE8866",  # ST65
  "#CC99BB",  # ST86
  "#AAAA00",  # ST29
  "#332288",  # ST258
  "#77AADD"   # ST512
), top_sts)

hh_cols <- c(
  "Liver abscess"            = "#994455",
  "Liver transplant"         = "#4477AA",
  "Other liver disease"      = "#228833",
  "Non-liver or unspecified" = "#887799",
  "Unknown"                  = "#DDDDDD")
kd$HH <- factor(kd$HH, levels=names(hh_cols))

row_ann <- rowAnnotation(
  ST = split_vec,
  `Host Health` = kd$HH,
  col = list(ST=st_cols, `Host Health`=hh_cols),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(6, "mm"),
  gap = unit(2, "mm"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

n_rows <- nrow(mat)
fig_height <- max(14, n_rows * 0.09 + 5)

# Use cell_fun to colour each cell by its locus colour
ht <- Heatmap(
  mat,
  name = "Virulence locus",
  col = c("0"="grey96", "1"="#5b5b5b"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(mat[i,j] == 1) {
      grid.rect(x, y, w, h,
                gp = gpar(fill=loci_colours[loci[j]], col="white", lwd=1))
    } else {
      grid.rect(x, y, w, h,
                gp = gpar(fill="grey96", col="white", lwd=1))
    }
  },
  top_annotation = top_ann,
  row_split = split_vec,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_side = "top",
  column_names_gp = gpar(fontsize=10, fontface="bold.italic"),
  column_names_rot = 30,
  column_title = NULL,
  row_title = st_labels[top_sts],
  row_title_gp = gpar(fontsize=10, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(4, "mm"),
  height = unit(n_rows * 2.5, "mm"),
  left_annotation = row_ann,
  right_annotation = right_ann,
  heatmap_legend_param = list(
    title = "Virulence locus",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey96","#5b5b5b")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)),
  border = TRUE,
  border_gp = gpar(col="grey85", lwd=0.3),
  use_raster = FALSE)

# Virulence score legend
vir_legend <- Legend(
  title = "Virulence score",
  col_fun = vir_col_fun,
  title_gp = gpar(fontsize=9, fontface="bold"),
  labels_gp = gpar(fontsize=8),
  at = c(0,1,2,3,4,5),
  direction = "horizontal")

# Loci colour legend
loci_legend <- Legend(
  title = "Virulence locus",
  labels = loci,
  legend_gp = gpar(fill=loci_colours[loci]),
  title_gp = gpar(fontsize=9, fontface="bold"),
  labels_gp = gpar(fontsize=8))

png(file.path(OUT,"Fig3_virulence_heatmap_v1.png"),
    width=12, height=fig_height, units="in", res=300)
draw(ht,
     merge_legend=FALSE,
     annotation_legend_list = list(loci_legend, vir_legend),
     padding=unit(c(8,12,8,15),"mm"),
     background="white")
dev.off()

pdf(file.path(OUT,"Fig3_virulence_heatmap_v1.pdf"), width=12, height=fig_height)
draw(ht,
     merge_legend=FALSE,
     annotation_legend_list = list(loci_legend, vir_legend),
     padding=unit(c(8,12,8,15),"mm"))
dev.off()

message(paste0("Fig3 v1 saved - ", n_rows, " isolates"))
