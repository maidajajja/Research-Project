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

# Categorical presence: 0=absent, 1=truncated/incomplete, 2=intact
# This is the key biological distinction
kleb$Yersiniabactin <- case_when(
  kleb[[ybt_col]] == "-"                     ~ 0,
  grepl("truncated|incomplete", kleb[[ybt_col]], ignore.case=TRUE) ~ 1,
  TRUE ~ 2)

kleb$Aerobactin <- case_when(
  kleb[[iuc_col]] == "-"                     ~ 0,
  grepl("incomplete", kleb[[iuc_col]], ignore.case=TRUE) ~ 1,
  TRUE ~ 2)

kleb$Salmochelin <- case_when(
  kleb[[iro_col]] == "-" ~ 0,
  TRUE ~ 2)

kleb$RmpADC <- case_when(
  kleb[[rmp_col]] == "-"                     ~ 0,
  grepl("truncated", kleb[[rmp_col]], ignore.case=TRUE) ~ 1,
  TRUE ~ 2)

kleb$rmpA2 <- case_when(
  kleb[[rmp2_col]] == "-"     ~ 0,
  grepl("\\*", kleb[[rmp2_col]]) ~ 1,  # truncated rmpA2 marked with *
  TRUE ~ 2)

kleb$Colibactin <- case_when(
  kleb[[clb_col]] == "-" ~ 0,
  TRUE ~ 2)

# hvKP status: score >=4
kleb$hvKP <- ifelse(kleb$vir_score >= 4, "hvKP (score>=4)", "Non-hvKP")

# Host Health
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

# 7 key STs ordered by decreasing virulence
top_sts <- c("ST23","ST11","ST65","ST86","ST29","ST258","ST512")
kd <- kleb %>% filter(ST %in% top_sts)
kd$ST_group <- factor(kd$ST, levels=top_sts)
kd <- kd %>% arrange(ST_group, desc(vir_score))
split_vec <- factor(as.character(kd$ST_group), levels=top_sts)

st_counts <- table(split_vec)
st_labels <- paste0(top_sts, " (n=", st_counts[top_sts], ")")
names(st_labels) <- top_sts

# Column order follows virulence score contribution
loci <- c("Yersiniabactin","Aerobactin","Salmochelin","RmpADC","rmpA2","Colibactin")
mat <- as.matrix(kd[, loci])
rownames(mat) <- kd$strain

cat("Distribution of values per locus:\n")
for(l in loci) cat(l, ":", table(mat[,l]), "\n")

# Prevalence of intact+truncated per locus
loci_prev <- colSums(mat > 0) / nrow(mat) * 100
loci_intact_prev <- colSums(mat == 2) / nrow(mat) * 100
cat("\nLoci prevalence (any):\n"); print(round(loci_prev,1))
cat("Loci prevalence (intact only):\n"); print(round(loci_intact_prev,1))

# Colour scheme: 3 states
# 0 = absent: grey96
# 1 = truncated/incomplete: light blue #88BBDD
# 2 = intact: dark blue #4477AA
cell_cols <- c("0"="grey96", "1"="#88BBDD", "2"="#4477AA")

# Prevalence bars - use intact prevalence, coloured dark blue
top_ann <- HeatmapAnnotation(
  `Prevalence\n(%)` = anno_barplot(
    cbind(loci_intact_prev, loci_prev - loci_intact_prev),
    gp = gpar(fill=c("#4477AA","#88BBDD"), col=NA),
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

# hvKP annotation on right
hvkp_cols <- c("hvKP (score>=4)"="#994455", "Non-hvKP"="grey90")
kd$hvKP <- factor(kd$hvKP, levels=names(hvkp_cols))

right_ann <- rowAnnotation(
  `hvKP\nstatus` = kd$hvKP,
  col = list(`hvKP\nstatus` = hvkp_cols),
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  annotation_name_side = "bottom",
  annotation_name_rot = 0,
  simple_anno_size = unit(6, "mm"),
  annotation_legend_param = list(
    title = "hvKP status",
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)))

# ST and Host Health
st_cols <- setNames(c(
  "#774411","#225522","#EE8866",
  "#CC99BB","#AAAA00","#332288","#77AADD"), top_sts)

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

ht <- Heatmap(
  mat,
  name = "Virulence locus",
  col = cell_cols,
  rect_gp = gpar(col="white", lwd=1.2),
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
    labels = c("Absent","Truncated/incomplete","Intact"),
    at = c(0,1,2),
    legend_gp = gpar(fill=c("grey96","#88BBDD","#4477AA")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)),
  border = TRUE,
  border_gp = gpar(col="grey85", lwd=0.3),
  use_raster = TRUE)

png(file.path(OUT,"Fig3_virulence_heatmap_v2.png"),
    width=12, height=fig_height, units="in", res=300)
draw(ht,
     merge_legend=TRUE,
     padding=unit(c(8,12,8,15),"mm"),
     background="white")
dev.off()

pdf(file.path(OUT,"Fig3_virulence_heatmap_v2.pdf"), width=12, height=fig_height)
draw(ht, merge_legend=TRUE, padding=unit(c(8,12,8,15),"mm"))
dev.off()

message(paste0("Fig3 v2 saved - ", n_rows, " isolates, 3-state: absent/truncated/intact"))
