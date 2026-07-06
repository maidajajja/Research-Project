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
amr$Accession <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", amr$Sample)
amr$Accession <- ifelse(grepl("^GC[AF]_", amr$Sample), amr$Accession, amr$Sample)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Accession_gca <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1",
                           trimws(meta[["Assembly Accession"]]))
meta$Accession_num <- trimws(as.character(meta[["Genome ID"]]))
meta$ST <- sub(".*\\.","", meta[["MLST"]])
meta$ST[is.na(meta$ST)|meta$ST==""] <- "Unknown"
meta$IsolationSource <- meta[["Isolation Source"]]
meta$IsolationSource[is.na(meta$IsolationSource)|meta$IsolationSource==""] <- "Unknown"
meta$Country <- meta[["Isolation Country"]]
meta$Country[is.na(meta$Country)|meta$Country==""|meta$Country=="NA"] <- "Unknown"

meta_lookup <- bind_rows(
  meta[, c("Accession_gca","ST","IsolationSource","Country")] %>%
    rename(Accession=Accession_gca),
  meta[, c("Accession_num","ST","IsolationSource","Country")] %>%
    rename(Accession=Accession_num)
) %>% filter(!is.na(Accession) & Accession != "") %>%
  distinct(Accession, .keep_all=TRUE)

amr2 <- left_join(amr, meta_lookup, by="Accession")
amr2 <- amr2[!is.na(amr2$Type) & amr2$Type == "AMR", ]
amr2$ST[is.na(amr2$ST)] <- "Unknown"

amr2$Source_clean <- case_when(
  grepl("liver abscess|hepatic abscess", amr2$IsolationSource, ignore.case=TRUE) ~ "Liver abscess",
  grepl("blood|bacteremia|bacteraemia", amr2$IsolationSource, ignore.case=TRUE) ~ "Blood",
  grepl("drain|pus", amr2$IsolationSource, ignore.case=TRUE) ~ "Drainage/Pus",
  grepl("bile|biliary|gallbladder", amr2$IsolationSource, ignore.case=TRUE) ~ "Bile/Biliary",
  grepl("liver", amr2$IsolationSource, ignore.case=TRUE) ~ "Liver",
  amr2$IsolationSource == "Unknown" ~ "Unknown",
  TRUE ~ "Other"
)

top_genes <- names(sort(table(amr2$Element.symbol), decreasing=TRUE))[1:25]
amr_filt <- amr2[amr2$Element.symbol %in% top_genes, ]

pa <- amr_filt %>%
  distinct(Accession, Element.symbol) %>%
  mutate(present=1) %>%
  pivot_wider(names_from=Element.symbol, values_from=present, values_fill=0)

sample_meta <- amr_filt %>%
  distinct(Accession, ST, Source_clean, Country) %>%
  group_by(Accession) %>% slice(1) %>% ungroup()

pa2 <- left_join(pa, sample_meta, by="Accession")
pa2$ST[is.na(pa2$ST)] <- "Unknown"
pa2$Source_clean[is.na(pa2$Source_clean)] <- "Unknown"
pa2$Country[is.na(pa2$Country)] <- "Unknown"

# ST ordering by ASCENDING mean gene count (fewest genes = top, most = bottom)
st_mean_genes <- pa2 %>%
  mutate(total=rowSums(across(all_of(top_genes)))) %>%
  group_by(ST) %>%
  summarise(mean_genes=mean(total), .groups="drop") %>%
  filter(ST != "Unknown") %>%
  arrange(mean_genes)

top_sts <- head(st_mean_genes$ST, 8)
st_levels <- c(paste0("ST", top_sts), "Other")
pa2$ST_group <- factor(
  ifelse(pa2$ST %in% top_sts, paste0("ST", pa2$ST), "Other"),
  levels=st_levels)

# within each ST: ascending gene count = stepped pattern
pa2 <- pa2 %>% arrange(ST_group, rowSums(across(all_of(top_genes))))

mat <- as.matrix(pa2[, top_genes])
rownames(mat) <- pa2$Accession
gene_order <- names(sort(colSums(mat), decreasing=TRUE))
mat <- mat[, gene_order]

cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#4D4D4D","#44AA99","grey75")
st_cols <- setNames(cbf_pal[1:length(st_levels)], st_levels)
src_cols <- c("Liver abscess"="#D55E00","Blood"="#CC79A7",
              "Drainage/Pus"="#56B4E9","Bile/Biliary"="#E69F00",
              "Liver"="#009E73","Other"="#4D4D4D","Unknown"="grey85")
cty_top <- names(sort(table(pa2$Country[pa2$Country!="Unknown"]),
                      decreasing=TRUE))[1:6]
cty_cols <- setNames(c("#0072B2","#E69F00","#009E73",
                       "#D55E00","#CC79A7","#56B4E9","grey75"),
                     c(cty_top,"Other"))
pa2$Country_group <- ifelse(pa2$Country %in% cty_top, pa2$Country, "Other")

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
vir_col <- grep("virulence_score__virulence_score", colnames(kleb), value=TRUE)[1]
kleb$res_score <- suppressWarnings(as.numeric(kleb[[res_col]]))
kleb$vir_score <- suppressWarnings(as.numeric(kleb[[vir_col]]))
score_map <- kleb[, c("strain","res_score","vir_score")]
pa2 <- left_join(pa2, score_map, by=c("Accession"="strain"))
pa2$res_score[is.na(pa2$res_score)] <- 0
pa2$vir_score[is.na(pa2$vir_score)] <- 0

row_ann <- rowAnnotation(
  ST = pa2$ST_group,
  `Vir. Score` = pa2$vir_score,
  `Res. Score` = pa2$res_score,
  `Isolation Source` = pa2$Source_clean,
  col = list(
    ST = st_cols,
    `Vir. Score` = colorRamp2(c(0,2,4,5), c("grey95","#F0E442","#E69F00","#D55E00")),
    `Res. Score` = colorRamp2(c(0,1,2,3), c("grey95","#56B4E9","#0072B2","#000033")),
    `Isolation Source` = src_cols
  ),
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(5, "mm"),
  gap = unit(1, "mm")
)

ht <- Heatmap(
  mat,
  name = "AMR gene",
  col = c("0"="grey95", "1"="black"),
  row_split = pa2$ST_group,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_distance_columns = "binary",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize=8, fontface="italic"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize=10, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),
  left_annotation = row_ann,
  heatmap_legend_param = list(
    title = "AMR gene",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey95","black")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)
  ),
  border = TRUE,
  border_gp = gpar(col="grey70", lwd=0.5),
  use_raster = FALSE
)

png(file.path(OUT,"Fig2_AMR_heatmap.png"),
    width=14, height=12, units="in", res=600)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,10),"mm"),
     background="white")
dev.off()

pdf(file.path(OUT,"Fig2_AMR_heatmap.pdf"), width=14, height=12)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,10),"mm"))
dev.off()

message("Fig2 saved at 600 DPI")
