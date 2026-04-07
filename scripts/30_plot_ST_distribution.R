library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Read Kleborate results for ST data
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

# Get ST column
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
kleb$ST <- kleb[[st_col]]

# Clean ST names - remove LV notation
kleb$ST <- gsub("-.*LV", "", kleb$ST)

# Count STs
st_counts <- kleb %>%
  group_by(ST) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Highlight major STs
major_sts <- c("ST23", "ST86", "ST11", "ST258", "ST65", "ST29", "ST512", "ST76", "ST60", "ST380")
st_counts$colour <- ifelse(st_counts$ST %in% major_sts, st_counts$ST, "Other")

# Custom colours
colours <- c(
  "ST23" = "#00BCD4", "ST86" = "#1565C0", "ST11" = "#FF6F00",
  "ST258" = "#E65100", "ST65" = "#6A1B9A", "ST29" = "#AD1457",
  "ST512" = "#2E7D32", "ST76" = "#00838F", "ST60" = "#558B2F",
  "ST380" = "#795548", "Other" = "#BDBDBD"
)

p <- ggplot(st_counts, aes(x = reorder(ST, -count), y = count, fill = colour)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = colours, name = "Sequence Type") +
  labs(title = "Sequence Type Distribution of 234 K. pneumoniae Isolates",
       x = "Sequence Type",
       y = "Number of Isolates") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank())

ggsave("/scratch/users/k22017808/KP_Research_Project/plots/ST_distribution.png",
       plot = p, width = 14, height = 6, dpi = 300)
cat("ST distribution plot saved\n")
