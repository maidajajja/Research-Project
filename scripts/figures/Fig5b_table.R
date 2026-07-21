suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"

mlst <- read.table("/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST/mlst_results.tsv",
                   sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(mlst) <- c("File","Scheme","ST","a1","a2","a3","a4","a5","a6","a7")
mlst$Sample <- gsub(".fasta$|.fna$|.fa$","", mlst$File)
mlst$ST <- as.character(mlst$ST)
mlst$ST[mlst$ST=="-"] <- NA

meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample  <- as.character(meta[["Genome ID"]])
meta$Country <- meta[["Isolation Country"]]
meta <- meta[, c("Sample","Country")]

df <- merge(mlst[,c("Sample","ST")], meta, by="Sample", all.x=TRUE)
df <- df[!is.na(df$ST) & !is.na(df$Country) & 
         df$Country != "" & df$Country != "NA" & df$Country != "Unknown", ]

top_sts <- names(sort(table(df$ST), decreasing=TRUE))[1:10]
st_labels <- paste0("ST", top_sts)

df_cty <- df %>%
  filter(ST %in% top_sts) %>%
  mutate(ST_label = paste0("ST", ST))

top_countries <- names(sort(table(df_cty$Country), decreasing=TRUE))[1:8]
df_cty <- df_cty %>% filter(Country %in% top_countries)

tbl <- df_cty %>%
  group_by(Country, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  pivot_wider(names_from=ST_label, values_from=n, values_fill=0)

# ensure all ST columns present
for(st in st_labels){ if(!st %in% colnames(tbl)) tbl[[st]] <- 0 }
tbl <- tbl[, c("Country", st_labels)]
tbl$Total <- rowSums(tbl[, st_labels])
tbl <- tbl %>% arrange(desc(Total))

# Convert to long format for ggplot table
col_names <- c("Country", st_labels, "Total")
tbl_display <- tbl[, col_names]

# Build plot-table using geom_tile + geom_text
n_rows <- nrow(tbl_display)
n_cols <- length(col_names)

cell_df <- expand.grid(row=1:n_rows, col=1:n_cols)
cell_df$value <- ""
for(i in 1:n_rows){
  for(j in 1:n_cols){
    v <- tbl_display[i, j]
    cell_df$value[cell_df$row==i & cell_df$col==j] <- as.character(v)
  }
}

# header row
header_df <- data.frame(row=0, col=1:n_cols, value=col_names)

all_df <- rbind(header_df, cell_df)
all_df$is_header <- all_df$row == 0
all_df$is_total_col <- all_df$col == n_cols
all_df$is_country_col <- all_df$col == 1
all_df$fill <- ifelse(all_df$is_header, "grey20",
               ifelse(all_df$row %% 2 == 0, "grey95", "white"))

figB <- ggplot(all_df, aes(x=col, y=-row)) +
  geom_tile(aes(fill=fill), colour="white", linewidth=0.5) +
  geom_text(aes(label=value,
                fontface=ifelse(is_header | is_total_col, "bold", "plain"),
                colour=ifelse(is_header, "white", "grey15")),
            size=3.5) +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(expand=c(0,0.5)) +
  scale_y_continuous(expand=c(0,0.5)) +
  theme_void() +
  theme(plot.margin=margin(10,10,10,10))

ggsave(file.path(OUT,"Fig5b_ST_country_table.png"), figB, 
       width=14, height=3, dpi=300, bg="white")
message("Table saved")
