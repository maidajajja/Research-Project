suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
dir.create(OUT, showWarnings=FALSE)

# Load MOBsuite data
mob_dir <- "/scratch/users/k22017808/KP_Research_Project/11_Plasmids/MOBsuite"
samples <- list.dirs(mob_dir, recursive=FALSE, full.names=FALSE)
samples <- samples[samples != "__tmp"]

all_mob <- lapply(samples, function(s) {
  f <- file.path(mob_dir, s, "contig_report.txt")
  if(!file.exists(f)) return(NULL)
  tryCatch({
    d <- read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE,
                    fill=TRUE, quote="")
    d$sample_id <- s; d
  }, error=function(e) NULL)
})
all_mob <- bind_rows(Filter(Negate(is.null), all_mob))

plasmids <- all_mob %>%
  filter(molecule_type == "plasmid") %>%
  select(sample_id, rep_type.s.)

# Correct MOBsuite replicon names
key_replicons <- c("IncHI1B","IncFII","IncFIA","IncR","IncC")

# Split comma-separated values and check presence
for(rep in key_replicons) {
  plasmids[[rep]] <- sapply(plasmids$rep_type.s., function(x) {
    if(is.na(x)) return(FALSE)
    any(trimws(unlist(strsplit(x, ","))) == rep)
  })
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

# Replicon prevalence by ST
st_rep <- df %>%
  filter(ST %in% top_sts) %>%
  group_by(ST) %>%
  summarise(
    n = n(),
    across(all_of(key_replicons), ~round(sum(.)/n()*100,1)),
    .groups="drop") %>%
  pivot_longer(cols=all_of(key_replicons),
               names_to="Replicon", values_to="Prevalence")

cat("ST counts in MOBsuite:\n")
print(df %>% filter(ST %in% top_sts) %>% count(ST))
cat("\nReplicon prevalence summary:\n")
print(st_rep %>% filter(Prevalence > 0) %>% arrange(Replicon, desc(Prevalence)))

st_rep$ST_label <- factor(paste0("ST", st_rep$ST),
                           levels=rev(paste0("ST", top_sts)))
st_rep$Replicon <- factor(st_rep$Replicon, levels=key_replicons)

# ST n labels
st_n <- df %>% filter(ST %in% top_sts) %>%
  group_by(ST) %>% summarise(n=n(), .groups="drop")
st_labels_n <- setNames(paste0("ST", st_n$ST, " (n=", st_n$n, ")"),
                        paste0("ST", st_n$ST))

# Paul Tol muted palette
rep_cols <- c(
  "IncHI1B" = "#4477AA",  # medium blue - virulence plasmid
  "IncFII"  = "#CC6677",  # muted rose - resistance plasmid
  "IncFIA"  = "#EEBB88",  # muted orange
  "IncR"    = "#BBBBBB",  # grey
  "IncC"    = "#997700"   # dark gold - carbapenemase associated
)

# ── Fig4a: Replicon prevalence by ST ─────────────────────────────────────────
fig4a <- ggplot(st_rep %>% filter(Prevalence > 0),
                aes(x=Prevalence, y=ST_label, fill=Replicon)) +
  geom_col(position="dodge", width=0.75, colour="white", linewidth=0.2) +
  geom_text(aes(label=paste0(Prevalence,"%")),
            position=position_dodge(width=0.75),
            hjust=-0.1, size=3, colour="grey30") +
  scale_fill_manual(values=rep_cols, name="Replicon type",
                    drop=FALSE) +
  scale_x_continuous(expand=expansion(mult=c(0,0.18)),
                     limits=c(0,115),
                     breaks=seq(0,100,25),
                     labels=function(x) paste0(x,"%")) +
  scale_y_discrete(labels=st_labels_n) +
  labs(x="Prevalence within ST (%)",
       y="Sequence Type",
       caption="Replicon typing by MOBsuite v3.1.9. Isolates may carry multiple replicon types.\nIncHI1B = virulence-associated conjugative plasmid; IncFII = resistance-associated plasmid.") +
  theme_classic(base_size=12) +
  theme(
    axis.title = element_text(face="bold", size=11),
    axis.text  = element_text(size=10, colour="grey20"),
    panel.grid.major.x = element_line(colour="grey90", linewidth=0.4),
    legend.title = element_text(face="bold", size=9),
    legend.text  = element_text(size=8.5),
    legend.position = "right",
    legend.background = element_blank(),
    plot.caption = element_text(size=7, colour="grey40", hjust=0),
    plot.margin = margin(10,20,10,10))

ggsave(file.path(OUT,"Fig4a_replicon_by_ST_v4.png"),
       fig4a, width=12, height=7, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4a_replicon_by_ST_v4.pdf"),
       fig4a, width=12, height=7)

# ── Fig4b: IncHI1B in ST11 convergent vs non-convergent ──────────────────────
st11 <- df %>%
  filter(ST == "11") %>%
  mutate(
    Group = ifelse(Convergent,
                   "Convergent\n(Vir>=4, Res>=2)",
                   "Non-convergent"),
    IncHI1B_status = ifelse(IncHI1B, "Present", "Absent"))

st11_summary <- st11 %>%
  group_by(Group, IncHI1B_status) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(Group) %>%
  mutate(pct=round(n/sum(n)*100,1), total=sum(n))

tbl <- table(st11$Convergent, st11$IncHI1B)
cat("\nST11 IncHI1B vs Convergent table:\n"); print(tbl)
ft <- fisher.test(tbl)
cat("Fisher's exact p =", round(ft$p.value,3), "\n")
pval_label <- paste0("Fisher's exact: p = ", round(ft$p.value,3),
                     ifelse(ft$p.value<0.05," *"," (n.s.)"))

group_n <- st11 %>% group_by(Group) %>% summarise(n=n())
group_labels <- setNames(paste0(group_n$Group,"\n(n=",group_n$n,")"),
                         group_n$Group)
st11_summary$Group_label <- group_labels[st11_summary$Group]
st11_summary$IncHI1B_status <- factor(st11_summary$IncHI1B_status,
                                       levels=c("Absent","Present"))

fig4b <- ggplot(st11_summary,
                aes(x=Group_label, y=pct, fill=IncHI1B_status)) +
  geom_col(width=0.5, colour="white", linewidth=0.4) +
  geom_text(aes(label=paste0(n,"\n(",pct,"%)"),
                colour=IncHI1B_status),
            position=position_stack(vjust=0.5),
            size=3.8, fontface="bold") +
  annotate("text", x=1.5, y=108,
           label=pval_label,
           size=3.5, colour="grey30", fontface="italic") +
  scale_fill_manual(values=c("Absent"="grey88","Present"="#4477AA"),
                    name="IncHI1B plasmid") +
  scale_colour_manual(values=c("Absent"="grey30","Present"="white"),
                      guide="none") +
  scale_y_continuous(limits=c(0,115),
                     breaks=seq(0,100,25),
                     labels=function(x) paste0(x,"%"),
                     expand=expansion(mult=c(0,0.05))) +
  labs(x=NULL,
       y="Proportion of ST11 isolates (%)",
       title="IncHI1B plasmid carriage in ST11",
       subtitle="Convergent vs non-convergent isolates",
       caption="Fisher's exact test. Low power due to small non-convergent group (n=4).\nConvergence defined as virulence score >=4 AND resistance score >=2 (Kleborate v3.2.4).") +
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

ggsave(file.path(OUT,"Fig4b_IncHI1B_ST11_v4.png"),
       fig4b, width=7, height=6, dpi=300, bg="white")
ggsave(file.path(OUT,"Fig4b_IncHI1B_ST11_v4.pdf"),
       fig4b, width=7, height=6)

message("Fig4 v4 saved")
