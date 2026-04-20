suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# ── Load data ────────────────────────────────────────────────────────────────
mlst <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST/mlst_results.tsv",
                   sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(mlst) <- c("File","Scheme","ST","a1","a2","a3","a4","a5","a6","a7")
mlst$Sample <- gsub(".fasta$|.fna$|.fa$","", mlst$File)
mlst$ST <- as.character(mlst$ST)
mlst$ST[mlst$ST=="-"] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Year    <- as.numeric(meta[["Collection Year"]])
meta$Country <- meta[["Isolation Country"]]
meta <- meta[, c("Sample","Year","Country")]

df <- merge(mlst[,c("Sample","ST")], meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$ST) & !is.na(df$Year) & df$Year > 2000, ]
df <- df[!is.na(df$Country) & df$Country != "" & 
         df$Country != "NA" & df$Country != "Unknown", ]

top_sts <- names(sort(table(df$ST), decreasing=TRUE))[1:10]
st_labels <- paste0("ST", top_sts)

# CBF palette (Wong 2011)
cbf_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
             "#D55E00","#CC79A7","#4D4D4D","#888888","#44AA99")
st_cols <- setNames(cbf_pal[1:length(st_labels)], st_labels)

pub_theme <- theme_classic(base_size=14) +
  theme(axis.title=element_text(face="bold", size=13),
        axis.text=element_text(size=11, colour="grey20"),
        legend.title=element_text(face="bold", size=11),
        legend.text=element_text(size=10),
        plot.margin=margin(12,24,12,12))

# ── Fig5a: Bubble plot ───────────────────────────────────────────────────────
st_year <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = factor(paste0("ST", ST), levels=rev(st_labels))) %>%
  group_by(Year, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  filter(n > 0)

actual_years <- sort(unique(st_year$Year))

# text colour: white for dark bubbles, dark for yellow (ST198/F0E442)
st_year$text_col <- ifelse(st_year$ST_label == "ST198", "grey20", "white")

figA <- ggplot(st_year, aes(x=Year, y=ST_label, size=n, fill=ST_label)) +
  geom_point(shape=21, colour="grey40", stroke=0.5, alpha=0.92) +
  geom_text(aes(label=n, colour=text_col), size=3.2, fontface="bold") +
  scale_size_continuous(range=c(7, 20), name="No. isolates",
                        breaks=c(1,3,6,9,12)) +
  scale_fill_manual(values=st_cols, guide="none") +
  scale_colour_identity() +
  scale_x_continuous(breaks=actual_years, labels=actual_years,
                     expand=expansion(add=c(0.3, 0.8))) +
  labs(x="Collection Year", y="Sequence Type") +
  pub_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_line(colour="grey85", linewidth=0.5),
        panel.grid.minor=element_blank(),
        legend.key=element_rect(fill="grey90", colour=NA))

ggsave(file.path(OUT,"Fig5a_ST_over_time.png"), figA, 
       width=9, height=6.5, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time.pdf"), figA, 
       width=9, height=6.5)

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

# relative column widths: Country=2, ST cols=1 each, Total=1.2
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
# replace zeros with em dash except Country col
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

# tile widths per column
all_df$tile_w <- col_widths[all_df$col]

figB <- ggplot(all_df) +
  geom_tile(aes(x=x, y=y, width=tile_w, height=0.85, fill=fill),
            colour="white", linewidth=0.4) +
  # outer border
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

message("Fig5 final saved")
