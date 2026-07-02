suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# NOTE: switched from the standalone mlst tool output to Kleborate's ST
# column - see Fig5_epidemiology_v7_klebfix.R for full rationale. The
# standalone mlst tool only covered 82/229 genomes.
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","", kleb[[st_col]])
kleb$ST <- sub("^ST", "", kleb$ST)
kleb$ST[kleb$ST==""] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)

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

meta$Year    <- as.numeric(meta[["Collection Year"]])
meta$Country <- meta[["Isolation Country"]]
meta <- meta[!is.na(meta$Sample), c("Sample","Year","Country")]

st_map <- data.frame(Sample=kleb$strain, ST=kleb$ST, stringsAsFactors=FALSE)

df <- merge(st_map, meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$ST) & !is.na(df$Year) & df$Year > 2000, ]
df <- df[!is.na(df$Country) & df$Country != "" &
         df$Country != "NA" & df$Country != "Unknown", ]

top_sts <- c("23","11","86","258","65","29","512")
verified_n <- c("23"=76, "11"=21, "86"=22, "258"=9, "65"=9, "29"=8, "512"=6)
st_labels <- paste0("ST", top_sts)
st_labels_n <- setNames(paste0("ST", names(verified_n), " (n=", verified_n, ")"), paste0("ST", names(verified_n)))

cbf_pal <- c("#774411","#225522","#CC99BB","#332288","#EE8866","#AAAA00","#77AADD",
             "#D55E00","#CC79A7","#4D4D4D","#888888","#44AA99")
st_cols <- setNames(cbf_pal[1:length(st_labels)], st_labels)

pub_theme <- theme_classic(base_size=14) +
  theme(axis.title=element_text(face="bold", size=13),
        axis.text=element_text(size=11, colour="grey20"),
        legend.title=element_text(face="bold", size=11),
        legend.text=element_text(size=10),
        plot.margin=margin(12,24,12,12))

# ── Fig5a: Bubble plot (superseded by Fig5a_fix3.R, kept for table panel only) ──
st_year <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = factor(paste0("ST", ST), levels=rev(st_labels))) %>%
  group_by(Year, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  filter(n > 0)

# ── Fig5b: Country table ─────────────────────────────────────────────────────
df_cty <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = paste0("ST", ST))

top_countries <- names(sort(table(df_cty$Country), decreasing=TRUE))[1:8]
df_cty <- df_cty %>% filter(Country %in% top_countries)

tbl <- df_cty %>%
  group_by(Country, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  pivot_wider(names_from=ST_label, values_from=n, values_fill=0)

for(st in st_labels){ if(!st %in% colnames(tbl)) tbl[[st]] <- 0 }
tbl <- tbl[, c("Country", st_labels)]
tbl$Total <- rowSums(tbl[, st_labels])
tbl <- tbl %>% arrange(desc(Total))

col_names <- c("Country", st_labels, "Total")
tbl_display <- tbl[, col_names]

n_rows <- nrow(tbl_display)
n_cols <- length(col_names)

col_widths <- c(2, rep(1, length(st_labels)), 1.2)
col_positions <- cumsum(c(0, col_widths))
col_centres <- col_positions[-length(col_positions)] + col_widths/2
total_width <- sum(col_widths)

cell_df <- expand.grid(row=seq_len(n_rows), col=seq_len(n_cols))
cell_df$value <- ""
cell_df$raw_n <- 0
for(i in seq_len(n_rows)){
  for(j in seq_len(n_cols)){
    v <- tbl_display[i, j]
    cell_df$value[cell_df$row==i & cell_df$col==j] <- as.character(v)
    cell_df$raw_n[cell_df$row==i & cell_df$col==j] <- suppressWarnings(as.numeric(v))
  }
}
cell_df$display <- ifelse(
  cell_df$col > 1 & !is.na(cell_df$raw_n) & cell_df$raw_n == 0,
  "\u2014", cell_df$value)

header_df <- data.frame(row=0, col=seq_len(n_cols),
                        value=col_names, display=col_names,
                        raw_n=NA)

all_df <- rbind(header_df, cell_df)
all_df$x <- col_centres[all_df$col]
all_df$y <- -all_df$row
all_df$is_header <- all_df$row == 0
all_df$is_total_col <- all_df$col == n_cols
all_df$is_country_col <- all_df$col == 1
all_df$is_dash <- all_df$display == "\u2014"

all_df$fill <- ifelse(all_df$is_header, "grey20",
               ifelse(all_df$row %% 2 == 0, "#F5F5F5", "white"))
all_df$text_colour <- ifelse(all_df$is_header, "white",
                      ifelse(all_df$is_dash, "grey70", "grey15"))
all_df$bold <- ifelse(all_df$is_header | all_df$is_total_col, "bold", "plain")

all_df$tile_w <- col_widths[all_df$col]

figB <- ggplot(all_df) +
  geom_tile(aes(x=x, y=y, width=tile_w, height=0.85, fill=fill),
            colour="white", linewidth=0.4) +
  annotate("rect", xmin=0, xmax=total_width,
           ymin=-(n_rows+0.425), ymax=0.425,
           fill=NA, colour="grey40", linewidth=0.6) +
  geom_text(aes(x=x, y=y, label=display,
                colour=text_colour, fontface=bold),
            size=3.4) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_void() +
  theme(plot.margin=margin(8,8,8,8))

ggsave(file.path(OUT,"Fig5b_ST_country_table.png"), figB,
       width=15, height=2.8, dpi=300, bg="white")

message("Fig5 final v2 saved - switched to Kleborate ST")
