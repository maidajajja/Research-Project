suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# ── Load MOBsuite data ────────────────────────────────────────────────────────
mob_dir <- "/scratch/users/k22017808/KP_Research_Project/11_Plasmids/MOBsuite"
samples <- list.dirs(mob_dir, recursive=FALSE, full.names=FALSE)
samples <- samples[samples != "__tmp"]

all_mob <- lapply(samples, function(s) {
  f <- file.path(mob_dir, s, "contig_report.txt")
  if(!file.exists(f)) return(NULL)
  tryCatch({
    d <- read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                    fill=TRUE, quote="")
    d$sample_id <- s
    d
  }, error=function(e) NULL)
})
all_mob <- bind_rows(Filter(Negate(is.null), all_mob))

# Get plasmid rows
plasmids <- all_mob %>%
  filter(molecule_type == "plasmid") %>%
  select(sample_id, primary_cluster_id, rep_type.s., predicted_mobility)

# Key replicon types to track
key_replicons <- c("IncHI1B","IncFIB(K)","ColRNAI","IncFII(K)",
                   "IncFII(pHN7A8)","IncR","IncFIB(pKPHS1)")

for(rep in key_replicons) {
  plasmids[[rep]] <- grepl(rep, plasmids$rep_type.s., ignore.case=TRUE)
}

# Per sample replicon presence
sample_replicons <- plasmids %>%
  group_by(sample_id) %>%
  summarise(across(all_of(key_replicons), any), .groups="drop")

# Load Kleborate
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt",
                   sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col  <- grep("mlst__ST$", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST  <- sub("^ST","", gsub("-.*","", kleb[[st_col]]))
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))
kleb$Convergent <- kleb$Vir >= 4 & kleb$Res >= 2

df <- left_join(sample_replicons,
                kleb[, c("strain","ST","Vir","Res","Convergent")],
                by=c("sample_id"="strain")) %>%
  filter(!is.na(ST))

top_sts <- c("23","11","86","258","65","29","512")

# ── Fig4a: Replicon prevalence by ST ─────────────────────────────────────────
st_rep <- df %>%
  filter(ST %in% top_sts) %>%
  group_by(ST) %>%
  summarise(
    n = n(),
    across(all_of(key_replicons), ~round(sum(.)/n()*100, 1)),
    .groups="drop") %>%
  pivot_longer(cols=all_of(key_replicons),
               names_to="Replicon", values_to="Prevalence")

st_rep$ST_label <- factor(paste0("ST", st_rep$ST),
                           levels=rev(paste0("ST", top_sts)))
st_rep$Replicon <- factor(st_rep$Replicon, levels=key_replicons)

# Add n to ST labels
st_n <- df %>% filter(ST %in% top_sts) %>%
  group_by(ST) %>% summarise(n=n(), .groups="drop")
st_labels_n <- paste0("ST", st_n$ST, " (n=", st_n$n, ")")
names(st_labels_n) <- paste0("ST", st_n$ST)

# Replicon colours - Paul Tol muted
rep_cols <- c(
  "IncHI1B"        = "#4477AA",
  "IncFIB(K)"      = "#44AA99",
  "ColRNAI"        = "#CCBB44",
  "IncFII(K)"      = "#EE6677",
  "IncFII(pHN7A8)" = "#AA3377",
  "IncR"           = "#BBBBBB",
  "IncFIB(pKPHS1)" = "#228833")

fig4a <- ggplot(st_rep, aes(x=Prevalence, y=ST_label, fill=Replicon)) +
  geom_col(position="dodge", width=0.7, colour="white", linewidth=0.2) +
  geom_text(aes(label=ifelse(Prevalence>0, paste0(Prevalence,"%"), "")),
            position=position_dodge(width=0.7),
            hjust=-0.1, size=2.5, colour="grey30") +
  scale_fill_manual(values=rep_cols, name="Replicon type") +
  scale_x_continuous(expand=expansion(mult=c(0,0.15)),
                     limits=c(0,115),
                     breaks=seq(0,100,25),
                     labels=function(x) paste0(x,"%")) +
  scale_y_discrete(labels=st_labels_n) +
  labs(x="Prevalence within ST (%)",
       y="Sequence Type",
       caption="Replicon typing by MOBsuite v3.1.9. Isolates may carry multiple replicon types.") +
  theme_classic(base_size=12) +
  theme(
    axis.title = element_text(face="bold", size=11),
    axis.text  = element_text(size=10, colour="grey20"),
    panel.grid.major.x = element_line(colour="grey90", linewidth=0.4),
    legend.title = element_text(face="bold", size=9),
    legend.text  = element_text(size=8),
    legend.position = "right",
    legend.background = element_blank(),
    plot.caption = element_text(size=7, colour="grey40", hjust=0),
    plot.margin = margin(10,15,10,10))

ggsave(file.path(OUT,"Fig4a_replicon_by_ST_v3.png"),
       fig4a, width=12, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4a_replicon_by_ST_v3.pdf"),
       fig4a, width=12, height=7)

# ── Fig4b: IncHI1B in ST11 convergent vs non-convergent ──────────────────────
st11 <- df %>%
  filter(ST == "11") %>%
  mutate(
    Group = ifelse(Convergent, "Convergent\n(Vir>=4, Res>=2)",
                               "Non-convergent"),
    IncHI1B_status = ifelse(IncHI1B, "Present", "Absent"))

st11_summary <- st11 %>%
  group_by(Group, IncHI1B_status) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(Group) %>%
  mutate(pct=round(n/sum(n)*100,1),
         total=sum(n))

# Fisher's test
tbl <- table(st11$Convergent, st11$IncHI1B)
ft <- fisher.test(tbl)
cat("Fisher's exact p =", round(ft$p.value, 3), "\n")
pval_label <- ifelse(ft$p.value < 0.05,
                     paste0("p = ", round(ft$p.value,3), " *"),
                     paste0("p = ", round(ft$p.value,3), " (n.s.)"))

# Group labels with n
group_n <- st11 %>% group_by(Group) %>% summarise(n=n())
group_labels <- setNames(
  paste0(group_n$Group, "\n(n=", group_n$n, ")"),
  group_n$Group)

st11_summary$Group_label <- group_labels[st11_summary$Group]
st11_summary$IncHI1B_status <- factor(st11_summary$IncHI1B_status,
                                       levels=c("Absent","Present"))

fig4b <- ggplot(st11_summary,
                aes(x=Group_label, y=pct, fill=IncHI1B_status)) +
  geom_col(width=0.5, colour="white", linewidth=0.4) +
  geom_text(aes(label=paste0(n, "\n(", pct, "%)")),
            position=position_stack(vjust=0.5),
            size=3.5, fontface="bold",
            colour=ifelse(st11_summary$IncHI1B_status=="Absent", "grey30", "white")) +
  annotate("text", x=1.5, y=105,
           label=pval_label,
           size=3.5, colour="grey30", fontface="italic") +
  scale_fill_manual(values=c("Absent"="grey88","Present"="#4477AA"),
                    name="IncHI1B plasmid") +
  scale_y_continuous(limits=c(0,115),
                     breaks=seq(0,100,25),
                     labels=function(x) paste0(x,"%"),
                     expand=expansion(mult=c(0,0.05))) +
  labs(x=NULL,
       y="Proportion of ST11 isolates (%)",
       title="IncHI1B plasmid carriage in ST11",
       subtitle="Convergent vs non-convergent isolates",
       caption=paste0("Fisher's exact test. Convergence = Vir score >=4 AND Res score >=2 (Kleborate).")) +
  theme_classic(base_size=12) +
  theme(
    axis.title = element_text(face="bold", size=11),
    axis.text  = element_text(size=10, colour="grey20"),
    panel.grid.major.y = element_line(colour="grey90", linewidth=0.4),
    legend.title = element_text(face="bold", size=9),
    legend.text  = element_text(size=8),
    plot.title   = element_text(face="bold", size=12),
    plot.subtitle = element_text(size=10, colour="grey40"),
    plot.caption = element_text(size=7, colour="grey40", hjust=0),
    plot.margin  = margin(10,15,10,10))

ggsave(file.path(OUT,"Fig4b_IncHI1B_ST11_convergent_v3.png"),
       fig4b, width=7, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4b_IncHI1B_ST11_convergent_v3.pdf"),
       fig4b, width=7, height=6)

cat("\nST11 IncHI1B summary:\n")
print(st11_summary)
message("Fig4 v3 saved")
