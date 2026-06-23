suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

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

kleb$Yersiniabactin <- case_when(
  kleb[[ybt_col]] == "-" ~ 0,
  grepl("truncated|incomplete", kleb[[ybt_col]], ignore.case=TRUE) ~ 1,
  TRUE ~ 2)
kleb$Aerobactin <- case_when(
  kleb[[iuc_col]] == "-" ~ 0,
  grepl("incomplete", kleb[[iuc_col]], ignore.case=TRUE) ~ 1,
  TRUE ~ 2)
kleb$Salmochelin <- case_when(
  kleb[[iro_col]] == "-" ~ 0, TRUE ~ 2)
kleb$RmpADC <- case_when(
  kleb[[rmp_col]] == "-" ~ 0,
  grepl("truncated", kleb[[rmp_col]], ignore.case=TRUE) ~ 1,
  TRUE ~ 2)
kleb$rmpA2 <- case_when(
  kleb[[rmp2_col]] == "-" ~ 0,
  grepl("\\*", kleb[[rmp2_col]]) ~ 1,
  TRUE ~ 2)
kleb$Colibactin <- case_when(
  kleb[[clb_col]] == "-" ~ 0, TRUE ~ 2)

top_sts <- c("ST23","ST11","ST65","ST86","ST29","ST258","ST512")
loci <- c("Yersiniabactin","Aerobactin","Salmochelin","RmpADC","rmpA2","Colibactin")

kd <- kleb %>% filter(ST %in% top_sts)
kd$ST <- factor(kd$ST, levels=rev(top_sts))  # reverse for ggplot y-axis

# Calculate per-ST per-locus summary
summary_df <- kd %>%
  group_by(ST) %>%
  summarise(
    n = n(),
    across(all_of(loci), list(
      prev  = ~round(sum(. > 0)/n()*100, 1),
      intact = ~round(sum(. == 2)/n()*100, 1),
      trunc  = ~round(sum(. == 1)/n()*100, 1)
    )),
    .groups="drop")

# Pivot to long format
prev_long <- summary_df %>%
  select(ST, n, ends_with("_prev")) %>%
  pivot_longer(cols=ends_with("_prev"),
               names_to="Locus",
               values_to="Prevalence") %>%
  mutate(Locus=gsub("_prev","",Locus))

intact_long <- summary_df %>%
  select(ST, ends_with("_intact")) %>%
  pivot_longer(cols=ends_with("_intact"),
               names_to="Locus",
               values_to="PctIntact") %>%
  mutate(Locus=gsub("_intact","",Locus))

plot_df <- left_join(prev_long, intact_long, by=c("ST","Locus"))
plot_df$Locus <- factor(plot_df$Locus, levels=loci)
plot_df$ST <- factor(plot_df$ST, levels=rev(top_sts))

# ST sample sizes for labels
st_n <- kd %>% group_by(ST) %>% summarise(n=n(), .groups="drop")
st_labels_vec <- paste0(st_n$ST, "\n(n=", st_n$n, ")")
names(st_labels_vec) <- st_n$ST

cat("Summary:\n"); print(plot_df)

# Paul Tol muted palette for intact % gradient
p <- ggplot(plot_df, aes(x=Locus, y=ST)) +
  # Background circle showing total prevalence
  geom_point(aes(size=Prevalence), colour="grey80", fill="grey80",
             shape=21, stroke=0) +
  # Foreground circle showing intact proportion
  geom_point(aes(size=PctIntact, fill=PctIntact),
             shape=21, stroke=0.5, colour="white") +
  # Add text labels for prevalence
  geom_text(aes(label=ifelse(Prevalence>0, paste0(Prevalence,"%"), "")),
            size=2.5, fontface="bold", colour="grey20") +
  scale_size_continuous(
    name="Prevalence (%)\n[total circle]",
    range=c(1,14),
    limits=c(0,100),
    breaks=c(25,50,75,100)) +
  scale_fill_gradientn(
    name="% Intact\n[fill colour]",
    colours=c("grey96","#88BBDD","#4477AA","#332288"),
    values=c(0,0.1,0.5,1),
    limits=c(0,100),
    breaks=c(0,25,50,75,100)) +
  scale_x_discrete(position="top") +
  scale_y_discrete(labels=st_labels_vec) +
  theme_minimal(base_size=11) +
  theme(
    panel.grid.major = element_line(colour="grey90", linewidth=0.3),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold.italic", size=10,
                                hjust=0, angle=-30),
    axis.text.y = element_text(face="bold", size=9),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size=8, face="bold"),
    legend.text = element_text(size=8),
    plot.margin = margin(20,10,10,10)) +
  guides(
    size = guide_legend(override.aes=list(fill="grey70", colour="grey70")),
    fill = guide_colourbar(barwidth=0.8, barheight=6))

ggsave(file.path(OUT,"Fig3b_virulence_bubble.png"),
       p, width=8, height=5, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig3b_virulence_bubble.pdf"),
       p, width=8, height=5)

message("Fig3b bubble plot saved")
