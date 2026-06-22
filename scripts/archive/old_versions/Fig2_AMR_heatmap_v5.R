suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
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
# (Same fix validated for Fig3_virulence_heatmap.R.)
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

meta$Country <- meta[["Isolation Country"]]
meta$Source  <- meta[["Isolation Source"]]
meta$Country[is.na(meta$Country)|meta$Country==""] <- "Unknown"
meta$Source[is.na(meta$Source)|meta$Source==""]   <- "Unknown"
meta$Source <- dplyr::recode(meta$Source,
  "liver abscess"="Liver abscess","hepatic abscess"="Liver abscess",
  "liver abscess drainage fluid"="Drainage/Pus",
  "liver abscess puncture fluid"="Drainage/Pus",
  "Liver abscess puncture fluid"="Drainage/Pus",
  "abscess drainage"="Drainage/Pus","liver drainage"="Drainage/Pus",
  "liver"="Liver","blood"="Blood","Blood"="Blood")
meta <- meta[!is.na(meta$Sample) & meta$Sample!="", c("Sample","Country","Source")]

st_col  <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST  <- gsub("ST", "", gsub("-.*","", kleb[[st_col]]))
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
  if(is.null(df)||nrow(df)==0||!"Element symbol" %in% colnames(df)) return(NULL)
  data.frame(Sample=sample, Gene=df[["Element symbol"]])
})
amr_all <- do.call(rbind, Filter(Negate(is.null), amr_list))

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
ann$ST[is.na(ann$ST)]         <- "Unknown"
ann$Source[is.na(ann$Source)] <- "Unknown"
ann$Vir[is.na(ann$Vir)]       <- 0
ann$Res[is.na(ann$Res)]       <- 0

# ── ST ordering by mean AMR gene count per isolate (ascending) ────────────────
ann$gene_count <- rowSums(mat[match(ann$Sample, rownames(mat)), ])
st_mean <- ann %>%
  filter(ST != "Unknown") %>%
  group_by(ST) %>%
  summarise(mean_genes=mean(gene_count), .groups="drop") %>%
  arrange(mean_genes)

top_sts_by_count <- names(sort(table(ann$ST[ann$ST!="Unknown"]), decreasing=TRUE))[1:7]
st_mean_filtered <- st_mean %>% filter(ST %in% top_sts_by_count)
top_sts_asc <- st_mean_filtered$ST

st_fixed <- c("23"="#E69F00","86"="#56B4E9","11"="#009E73",
              "258"="#F0E442","65"="#0072B2","29"="#D55E00",
              "512"="#CC79A7","Other"="grey80","Unknown"="grey90")
st_cols <- st_fixed[c(top_sts_asc,"Other","Unknown")]
names(st_cols) <- c(top_sts_asc,"Other","Unknown")
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

row_ann <- rowAnnotation(
  ST           = ann$ST_group,
  `Vir. Score` = ann$Vir,
  `Res. Score` = ann$Res,
  Source       = ann$Src_clean,
  col = list(
    ST           = st_cols,
    `Vir. Score` = colorRamp2(c(0,2,5), c("#FFF7BC","#FE9929","#CC4C02")),
    `Res. Score` = colorRamp2(c(0,1,3), c("#EDF8FB","#66C2A4","#00441B")),
    Source       = src_cols
  ),
  annotation_name_gp = gpar(fontsize=9, fontface="bold"),
  annotation_width   = unit(c(0.8,0.4,0.4,0.6),"cm"),
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
  column_title     = NULL,
  cluster_rows     = FALSE,
  clustering_distance_columns = "binary",
  clustering_method_columns   = "ward.D2",
  rect_gp = gpar(col="white", lwd=0.3),
  heatmap_legend_param = list(
    title="Present", labels=c("No","Yes"),
    at=c(0,1), legend_height=unit(2,"cm")
  )
)

png(file.path(OUT,"Fig2_AMR_heatmap.png"),
    width=13, height=14, units="in", res=600, bg="white")
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,5),"mm"))
dev.off()
pdf(file.path(OUT,"Fig2_AMR_heatmap.pdf"), width=13, height=14)
draw(ht, merge_legend=TRUE)
dev.off()
message("Fig2 v5 saved at 600 DPI - Kleborate ID-join fix applied")
