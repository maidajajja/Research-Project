suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load Kleborate
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST$", colnames(kleb), value=TRUE)[1]
kleb$ST <- sub("ST","", sub("-.*","", kleb[[st_col]]))

# Load metadata - use Kleborate_ID as join key
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv",
                 stringsAsFactors=FALSE, check.names=FALSE)
meta_lookup <- meta[meta$Kleborate_ID != "" & !is.na(meta$Kleborate_ID),
                    c("Kleborate_ID","Collection Year")]
colnames(meta_lookup) <- c("strain","Year")

cat("Kleborate strains:", nrow(kleb), "\n")
cat("Meta rows with Kleborate_ID:", nrow(meta_lookup), "\n")
cat("Matching:", sum(kleb$strain %in% meta_lookup$strain), "\n")

df <- left_join(kleb[,c("strain","ST")], meta_lookup, by="strain")

# Verified counts
verified_n <- c("23"=76,"11"=21,"86"=22,"258"=9,"65"=9,"29"=8,"512"=6)
top_sts <- names(verified_n)

# ST colours consistent with other figures
st_cols_base <- c(
  "23"="#774411","11"="#225522","86"="#CC99BB",
  "258"="#332288","65"="#EE8866","29"="#AAAA00","512"="#77AADD"
)

# Labels with n=
st_labels_n <- setNames(
  paste0("ST", names(verified_n), " (n=", verified_n, ")"),
  paste0("ST", names(verified_n))
)

# Filter to top STs, valid years only
st_year <- df %>%
  filter(ST %in% top_sts,
         !is.na(Year), Year != "NA", Year != "",
         !is.na(as.integer(Year)) & as.integer(Year) > 2000) %>%
  mutate(ST_label = factor(st_labels_n[paste0("ST", ST)],
                           levels=rev(st_labels_n))) %>%
  group_by(Year, ST_label) %>%
  summarise(n=n(), .groups="drop") %>%
  filter(n > 0)

cat("ST-year combinations:", nrow(st_year), "\n")

# Dynamic year positions
year_vals <- sort(unique(st_year$Year))
year_pos <- setNames(seq_along(year_vals), as.character(year_vals))
st_year$x_pos <- year_pos[as.character(st_year$Year)]

# Colour mapping
st_cols <- setNames(
  st_cols_base[sub("ST(\\d+).*","\\1", levels(st_year$ST_label))],
  levels(st_year$ST_label)
)

figA <- ggplot(st_year, aes(x=x_pos, y=ST_label, size=n, fill=ST_label)) +
  geom_point(shape=21, colour="grey30", stroke=0.6, alpha=0.9) +
  geom_text(aes(label=n), colour="white", size=3.2, fontface="bold") +
  scale_size_continuous(range=c(6,22), name="No. isolates",
                        breaks=c(1,5,10,20,27)) +
  scale_fill_manual(values=st_cols, guide="none") +
  scale_x_continuous(breaks=year_pos,
                     labels=as.character(year_vals),
                     expand=expansion(add=c(0.5,1.0))) +
  labs(x="Collection Year", y="Sequence Type") +
  theme_classic(base_size=14) +
  theme(
    axis.title=element_text(face="bold", size=13),
    axis.text=element_text(size=10, colour="grey20"),
    axis.text.x=element_text(angle=45, hjust=1),
    axis.text.y=element_text(margin=margin(r=8)),
    panel.grid.major=element_line(colour="grey88", linewidth=0.4),
    panel.grid.minor=element_blank(),
    legend.title=element_text(face="bold", size=11),
    legend.text=element_text(size=10),
    legend.background=element_blank(),
    plot.margin=margin(12,24,12,20)
  )

ggsave(file.path(OUT,"Fig5a_ST_over_time_v4.png"), figA,
       width=14, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig5a_ST_over_time_v4.pdf"), figA,
       width=14, height=7)

message("Fig5a v4 saved - ", nrow(st_year), " ST-year combinations")
# DEBUG
cat("Year column sample:\n")
print(head(meta_lookup$Year, 20))
cat("Year class:", class(meta_lookup$Year), "\n")
cat("df Year sample:\n")
print(head(df$Year, 20))
