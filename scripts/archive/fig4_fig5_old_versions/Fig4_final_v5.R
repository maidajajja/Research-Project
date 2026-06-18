suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

plasmid <- read.table("/scratch/users/k22017808/KP_Research_Project/11_Plasmids/plasmidfinder_results.tsv",
                      sep="\t", header=TRUE, stringsAsFactors=FALSE)
plasmid <- plasmid[plasmid$Plasmid != "None", ]

# extract accession OR numeric ID
plasmid$Accession <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1", plasmid$Sample)
# for non-GCA samples, accession = sample itself
plasmid$Accession <- ifelse(grepl("^GC[AF]_", plasmid$Sample),
                            plasmid$Accession, plasmid$Sample)

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Accession_GCA <- sub("^(GC[AF]_[0-9]+\\.[0-9]+).*","\\1",
                           trimws(meta[["Assembly Accession"]]))
# numeric genome ID as character for SRA genomes
meta$Accession_num <- trimws(as.character(meta[["Genome ID"]]))
meta$ST <- sub(".*\\.","", meta[["MLST"]])
meta$ST[is.na(meta$ST)|meta$ST==""] <- "Unknown"
meta$Country <- meta[["Isolation Country"]]
meta$Country[is.na(meta$Country)|meta$Country==""|meta$Country=="NA"] <- "Unknown"

# build two lookup tables and combine
meta_gca <- meta[!is.na(meta$Accession_GCA) & meta$Accession_GCA != "",
                 c("Accession_GCA","ST","Country")]
colnames(meta_gca)[1] <- "Accession"

meta_num <- meta[!is.na(meta$Accession_num) & meta$Accession_num != "",
                 c("Accession_num","ST","Country")]
colnames(meta_num)[1] <- "Accession"
# format numeric to match plasmid sample format e.g. 573.78610
meta_num$Accession <- formatC(as.numeric(meta_num$Accession), 
                               format="f", digits=5)
# actually just use as.character with enough precision
meta_num$Accession <- as.character(as.numeric(meta[["Genome ID"]]))

meta_lookup <- bind_rows(meta_gca, meta_num) %>%
  distinct(Accession, .keep_all=TRUE)

df <- left_join(plasmid, meta_lookup, by="Accession") %>%
  distinct(Accession, Plasmid, .keep_all=TRUE)
df$ST[is.na(df$ST)] <- "Unknown"
df$Country[is.na(df$Country)] <- "Unknown"

message(sprintf("Unique sample-plasmid combinations: %d", nrow(df)))
message(sprintf("Rows with ST assigned: %d", sum(df$ST != "Unknown")))
message("Top STs:")
print(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE)[1:10])
message("Top countries:")
print(sort(table(df$Country[df$Country!="Unknown"]), decreasing=TRUE)[1:8])

# Top 8 STs
top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]), decreasing=TRUE))[1:8]
st_labels_full <- paste0("ST", top_sts)
df$ST_group <- ifelse(df$ST %in% top_sts, paste0("ST", df$ST), "Other")
df$ST_group <- factor(df$ST_group, levels=c(st_labels_full, "Other"))

top_plasmids <- names(sort(table(df$Plasmid), decreasing=TRUE))[1:15]
df_filt <- df[df$Plasmid %in% top_plasmids, ]

cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#4D4D4D","grey75")
st_cols <- setNames(cbf_pal[1:length(levels(df$ST_group))],
                    levels(df$ST_group))

# ── Fig4a ────────────────────────────────────────────────────────────────────
plasmid_st <- df_filt %>%
  group_by(ST_group, Plasmid) %>%
  summarise(n=n(), .groups="drop")

plasmid_order <- plasmid_st %>%
  group_by(Plasmid) %>% summarise(total=sum(n)) %>%
  arrange(total) %>% pull(Plasmid)
plasmid_st$Plasmid <- factor(plasmid_st$Plasmid, levels=plasmid_order)

figA <- ggplot(plasmid_st, aes(x=n, y=Plasmid, fill=ST_group)) +
  geom_col(width=0.7, colour="white", linewidth=0.3) +
  scale_fill_manual(values=st_cols, name="Sequence Type") +
  scale_x_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(x="Number of isolates", y="Plasmid replicon type") +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold", size=12),
        axis.text=element_text(size=10, colour="grey20"),
        axis.text.y=element_text(face="italic"),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.4),
        legend.title=element_text(face="bold", size=10),
        legend.text=element_text(size=9),
        legend.background=element_blank(),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(OUT,"Fig4a_plasmid_by_ST.png"), figA,
       width=10, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4a_plasmid_by_ST.pdf"), figA, width=10, height=7)

# ── Fig4b: Table ─────────────────────────────────────────────────────────────
top_countries <- names(sort(table(
  df_filt$Country[df_filt$Country!="Unknown"]), decreasing=TRUE))[1:6]
df_cty <- df_filt[df_filt$Country %in% top_countries, ]
top10 <- names(sort(table(df_filt$Plasmid), decreasing=TRUE))[1:10]

tbl <- df_cty %>%
  group_by(Country, Plasmid) %>%
  summarise(n=n(), .groups="drop") %>%
  pivot_wider(names_from=Plasmid, values_from=n, values_fill=0)
for(p in top10){ if(!p %in% colnames(tbl)) tbl[[p]] <- 0 }
tbl <- tbl[, c("Country", top10[top10 %in% colnames(tbl)])]
tbl$Total <- rowSums(tbl[,-1])
tbl <- tbl %>% arrange(desc(Total))

col_names <- colnames(tbl)
n_rows <- nrow(tbl); n_cols <- length(col_names)
col_widths <- c(1.6, rep(1.3, n_cols-2), 0.9)
col_positions <- cumsum(c(0, col_widths))
col_centres <- col_positions[-length(col_positions)] + col_widths/2
total_width <- sum(col_widths)

cell_df <- expand.grid(row=seq_len(n_rows), col=seq_len(n_cols))
cell_df$value <- ""; cell_df$raw_n <- NA
for(i in seq_len(n_rows)) for(j in seq_len(n_cols)){
  v <- tbl[i,j]
  cell_df$value[cell_df$row==i&cell_df$col==j] <- as.character(v)
  cell_df$raw_n[cell_df$row==i&cell_df$col==j] <- suppressWarnings(as.numeric(v))
}
cell_df$display <- ifelse(cell_df$col>1 & !is.na(cell_df$raw_n) &
                          cell_df$raw_n==0, "\u2014", cell_df$value)
header_df <- data.frame(row=0, col=seq_len(n_cols), value=col_names,
                        display=col_names, raw_n=NA)
all_df <- rbind(header_df, cell_df)
all_df$x <- col_centres[all_df$col]
all_df$y <- -all_df$row
all_df$is_header <- all_df$row==0
all_df$is_total  <- all_df$col==n_cols
all_df$is_dash   <- all_df$display=="\u2014"
all_df$fill <- ifelse(all_df$is_header,"grey20",
               ifelse(all_df$row%%2==0,"#F5F5F5","white"))
all_df$text_colour <- ifelse(all_df$is_header,"white",
                      ifelse(all_df$is_dash,"grey70","grey15"))
all_df$bold <- ifelse(all_df$is_header|all_df$is_total,"bold","plain")
all_df$tile_w <- col_widths[all_df$col]
all_df$text_size <- ifelse(all_df$is_header, 2.8, 3.0)

figB <- ggplot(all_df) +
  geom_tile(aes(x=x, y=y, width=tile_w, height=0.85, fill=fill),
            colour="white", linewidth=0.4) +
  annotate("rect", xmin=0, xmax=total_width,
           ymin=-(n_rows+0.425), ymax=0.425,
           fill=NA, colour="grey40", linewidth=0.6) +
  geom_text(aes(x=x, y=y, label=display, colour=text_colour,
                fontface=bold, size=text_size)) +
  scale_fill_identity() + scale_colour_identity() + scale_size_identity() +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  theme_void() + theme(plot.margin=margin(8,8,8,8))

ggsave(file.path(OUT,"Fig4b_plasmid_by_country.png"), figB,
       width=15, height=2.8, dpi=300, bg="white")

message("Fig4 v5 final saved")
