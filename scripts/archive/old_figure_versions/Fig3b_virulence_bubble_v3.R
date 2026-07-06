suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"

kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

st_col   <- grep("mlst__ST$", colnames(kleb), value=TRUE)[1]
ybt_col  <- "klebsiella__ybst__Yersiniabactin"
clb_col  <- grep("cbst__Colibactin$", colnames(kleb), value=TRUE)[1]
iuc_col  <- grep("abst__Aerobactin$", colnames(kleb), value=TRUE)[1]
iro_col  <- grep("smst__Salmochelin$", colnames(kleb), value=TRUE)[1]
rmp_col  <- grep("rmst__RmpADC$", colnames(kleb), value=TRUE)[1]
rmp2_col <- grep("rmpa2__rmpA2$", colnames(kleb), value=TRUE)[1]

kleb$ST <- kleb[[st_col]]

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
kd$ST <- factor(kd$ST, levels=rev(top_sts))

summary_df <- kd %>%
  group_by(ST) %>%
  summarise(
    n = n(),
    across(all_of(loci), list(
      prev   = ~round(sum(. > 0)/n()*100, 1),
      intact = ~round(sum(. == 2)/n()*100, 1),
      trunc  = ~round(sum(. == 1)/n()*100, 1)
    )),
    .groups="drop")

prev_long <- summary_df %>%
  select(ST, n, ends_with("_prev")) %>%
  pivot_longer(ends_with("_prev"), names_to="Locus", values_to="Prevalence") %>%
  mutate(Locus=gsub("_prev","",Locus))

intact_long <- summary_df %>%
  select(ST, ends_with("_intact")) %>%
  pivot_longer(ends_with("_intact"), names_to="Locus", values_to="PctIntact") %>%
  mutate(Locus=gsub("_intact","",Locus))

plot_df <- left_join(prev_long, intact_long, by=c("ST","Locus"))
plot_df$Locus <- factor(plot_df$Locus, levels=loci)
plot_df$ST <- factor(plot_df$ST, levels=rev(top_sts))

# Only plot non-zero prevalence
plot_df_nonzero <- plot_df %>% filter(Prevalence > 0)

# Text colour — white for dark circles, dark grey for light circles
plot_df_nonzero <- plot_df_nonzero %>%
  mutate(text_col = ifelse(PctIntact > 40, "white", "grey20"))

st_n <- kd %>% group_by(ST) %>% summarise(n=n(), .groups="drop")
st_labels_vec <- paste0(as.character(st_n$ST), " (n=", st_n$n, ")")
names(st_labels_vec) <- as.character(st_n$ST)

p <- ggplot(plot_df_nonzero, aes(x=Locus, y=ST)) +
  geom_point(aes(size=Prevalence, fill=PctIntact),
             shape=21, stroke=0.4, colour="white") +
  # Adaptive text colour based on fill darkness
  geom_text(data=plot_df_nonzero %>% filter(text_col=="white"),
            aes(label=paste0(Prevalence,"%")),
            size=2.8, fontface="bold", colour="white") +
  geom_text(data=plot_df_nonzero %>% filter(text_col=="grey20"),
            aes(label=paste0(Prevalence,"%")),
            size=2.8, fontface="bold", colour="grey25") +
  scale_size_continuous(
    name="Prevalence (%)",
    range=c(3,16),
    limits=c(0,100),
    breaks=c(25,50,75,100)) +
  scale_fill_gradientn(
    name="% Intact",
    colours=c("#CCDDEE","#88BBDD","#4477AA","#332288"),
    limits=c(0,100),
    breaks=c(0,25,50,75,100)) +
  scale_x_discrete(position="top") +
  scale_y_discrete(labels=st_labels_vec) +
  labs(caption="Circle size = prevalence within ST. Fill colour = proportion of carriers with intact (fully functional) locus.\nST86 Colibactin absent in all isolates. Truncated/incomplete loci may have reduced virulence function.") +
  theme_minimal(base_size=11) +
  theme(
    panel.grid.major = element_line(colour="grey88", linewidth=0.4),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold.italic", size=10,
                                hjust=0, angle=-35, colour="grey20"),
    axis.text.y = element_text(face="bold", size=10, colour="grey20"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.title = element_text(size=9, face="bold"),
    legend.text = element_text(size=8),
    plot.caption = element_text(size=7, colour="grey40", hjust=0),
    plot.margin = margin(25,10,10,10),
    panel.background = element_rect(fill="grey98", colour=NA)) +
  guides(
    size = guide_legend(
      override.aes=list(fill="#6699CC", colour="white"),
      order=1),
    fill = guide_colourbar(
      barwidth=0.8, barheight=5, order=2))

ggsave(file.path(OUT,"Fig3b_virulence_bubble_v3.png"),
       p, width=9, height=5.5, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig3b_virulence_bubble_v3.pdf"),
       p, width=9, height=5.5)

message("Fig3b v3 saved")
