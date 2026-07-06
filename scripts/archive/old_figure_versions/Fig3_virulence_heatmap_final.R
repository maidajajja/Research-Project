suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Accession_gca <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1",
                           trimws(meta[["Assembly Accession"]]))
meta$Accession_num <- trimws(as.character(meta[["Genome ID"]]))
meta$ST <- sub(".*\\.","", meta[["MLST"]])
meta$ST[is.na(meta$ST)|meta$ST==""] <- "Unknown"
meta$HostHealth <- meta[["Host Health"]]
meta$HostHealth[is.na(meta$HostHealth)|meta$HostHealth==""] <- "Unknown"
meta$IsolationSource <- meta[["Isolation Source"]]
meta$IsolationSource[is.na(meta$IsolationSource)|meta$IsolationSource==""] <- "Unknown"

vir_cols <- c("klebsiella__ybst__Yersiniabactin",
              "klebsiella__cbst__Colibactin",
              "klebsiella__abst__Aerobactin",
              "klebsiella__smst__Salmochelin",
              "klebsiella__rmst__RmpADC",
              "klebsiella__rmpa2__rmpA2")
vir_cols <- vir_cols[vir_cols %in% colnames(kleb)]

kleb_vir <- kleb[, c("strain", vir_cols)]
for(col in vir_cols){
  kleb_vir[[col]] <- as.integer(
    trimws(kleb_vir[[col]]) != "-" & !is.na(kleb_vir[[col]]) &
    trimws(kleb_vir[[col]]) != "")
}

# join via accession number — format as character matching
kleb_vir$strain_char <- trimws(as.character(kleb_vir$strain))
meta$Accession_num_char <- trimws(as.character(meta$Accession_num))

df <- left_join(kleb_vir,
  meta[, c("Accession_num_char","Accession_gca","ST","HostHealth","IsolationSource")],
  by=c("strain_char"="Accession_num_char"))

message(sprintf("ST assigned: %d / %d", sum(!is.na(df$ST) & df$ST!="Unknown"), nrow(df)))

df$ST[is.na(df$ST)] <- "Unknown"
df$HostHealth[is.na(df$HostHealth)] <- "Unknown"
df$IsolationSource[is.na(df$IsolationSource)] <- "Unknown"

df$Health_clean <- case_when(
  grepl("healthy|asymptomatic", df$HostHealth, ignore.case=TRUE) ~ "Healthy",
  grepl("liver|hepatic|cirrhosis", df$HostHealth, ignore.case=TRUE) ~ "Liver disease",
  grepl("ill|sick|infect|disease|cancer|abscess", df$HostHealth, ignore.case=TRUE) ~ "Diseased",
  TRUE ~ "Unknown"
)

df$Source_clean <- case_when(
  grepl("liver abscess|hepatic abscess", df$IsolationSource, ignore.case=TRUE) ~ "Liver abscess",
  grepl("blood|bacteremia|bacteraemia", df$IsolationSource, ignore.case=TRUE) ~ "Blood",
  grepl("drain|pus", df$IsolationSource, ignore.case=TRUE) ~ "Drainage/Pus",
  grepl("bile|biliary|gallbladder", df$IsolationSource, ignore.case=TRUE) ~ "Bile/Biliary",
  grepl("liver", df$IsolationSource, ignore.case=TRUE) ~ "Liver",
  df$IsolationSource == "Unknown" ~ "Unknown",
  TRUE ~ "Other"
)

# ST ordering by ascending mean virulence locus count
st_mean <- df %>%
  mutate(total=rowSums(across(all_of(vir_cols)))) %>%
  group_by(ST) %>%
  summarise(mean_vir=mean(total), .groups="drop") %>%
  filter(ST != "Unknown") %>%
  arrange(mean_vir)

top_sts <- head(st_mean$ST, 8)
st_levels <- c(paste0("ST", top_sts), "Other")
df$ST_group <- factor(
  ifelse(df$ST %in% top_sts, paste0("ST", df$ST), "Other"),
  levels=st_levels)

# ascending within ST
df <- df %>% arrange(ST_group, rowSums(across(all_of(vir_cols))))

mat <- as.matrix(df[, vir_cols])
rownames(mat) <- df$strain
colnames(mat) <- c("Yersiniabactin","Colibactin","Aerobactin",
                   "Salmochelin","RmpADC","rmpA2")

cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#4D4D4D","#44AA99","grey75")
st_cols <- setNames(cbf_pal[1:length(st_levels)], st_levels)
health_cols <- c("Healthy"="#009E73","Liver disease"="#E69F00",
                 "Diseased"="#D55E00","Unknown"="grey85")
src_cols <- c("Liver abscess"="#D55E00","Blood"="#CC79A7",
              "Drainage/Pus"="#56B4E9","Bile/Biliary"="#E69F00",
              "Liver"="#009E73","Other"="#4D4D4D","Unknown"="grey85")

kleb$vir_score <- suppressWarnings(
  as.numeric(kleb$klebsiella_pneumo_complex__virulence_score__virulence_score))
df <- left_join(df, kleb[, c("strain","vir_score")], by="strain")
df$vir_score[is.na(df$vir_score)] <- 0

row_ann <- rowAnnotation(
  ST = df$ST_group,
  `Vir. Score` = df$vir_score,
  `Host Health` = df$Health_clean,
  `Isolation Source` = df$Source_clean,
  col = list(
    ST = st_cols,
    `Vir. Score` = colorRamp2(c(0,2,4,5), c("grey95","#F0E442","#E69F00","#D55E00")),
    `Host Health` = health_cols,
    `Isolation Source` = src_cols
  ),
  annotation_name_gp = gpar(fontsize=8, fontface="bold"),
  annotation_name_side = "bottom",
  simple_anno_size = unit(5, "mm"),
  gap = unit(1, "mm")
)

ht <- Heatmap(
  mat,
  name = "Virulence locus",
  col = c("0"="grey95", "1"="#0072B2"),
  row_split = df$ST_group,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize=10, fontface="bold"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize=10, fontface="bold"),
  row_title_rot = 0,
  row_gap = unit(2, "mm"),
  left_annotation = row_ann,
  heatmap_legend_param = list(
    title = "Virulence locus",
    labels = c("Absent","Present"),
    at = c(0,1),
    legend_gp = gpar(fill=c("grey95","#0072B2")),
    title_gp = gpar(fontsize=9, fontface="bold"),
    labels_gp = gpar(fontsize=8)
  ),
  border = TRUE,
  border_gp = gpar(col="grey70", lwd=0.5),
  use_raster = FALSE
)

png(file.path(OUT,"Fig3_virulence_heatmap.png"),
    width=10, height=12, units="in", res=600)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,10),"mm"),
     background="white")
dev.off()

pdf(file.path(OUT,"Fig3_virulence_heatmap.pdf"), width=10, height=12)
draw(ht, merge_legend=TRUE, padding=unit(c(5,5,5,10),"mm"))
dev.off()

message("Fig3 saved at 600 DPI")
