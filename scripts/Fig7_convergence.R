suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
OUT <- "/scratch/users/k22017808/KP_Research_Project/plots"
meta <- read.csv("/scratch/users/k22017808/KP_Research_Project/genomes.csv", stringsAsFactors=FALSE, check.names=FALSE)
meta$Sample <- as.character(meta[["Genome ID"]])
meta$Source <- meta[["Isolation Source"]]
meta$Health <- meta[["Host Health"]]
meta$Source[is.na(meta$Source)|meta$Source==""] <- "Unknown"
meta$Health[is.na(meta$Health)|meta$Health==""] <- "Unknown"
meta <- meta[, c("Sample","Source","Health")]
kleb <- read.table("/scratch/users/k22017808/KP_Research_Project/09_Kleborate/klebsiella_pneumo_complex_output.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
st_col <- grep("mlst__ST", colnames(kleb), value=TRUE)[1]
vir_col <- "klebsiella_pneumo_complex__virulence_score__virulence_score"
res_col <- grep("resistance_score__resistance_score", colnames(kleb), value=TRUE)[1]
kleb$ST <- gsub("-.*","",kleb[[st_col]])
kleb$Vir <- suppressWarnings(as.numeric(kleb[[vir_col]]))
kleb$Res <- suppressWarnings(as.numeric(kleb[[res_col]]))
df <- merge(data.frame(Sample=kleb$strain,ST=kleb$ST,Vir=kleb$Vir,Res=kleb$Res),meta,by="Sample",all.x=TRUE)
df$Source[is.na(df$Source)] <- "Unknown"
df$ST[is.na(df$ST)] <- "Unknown"
df <- df[!is.na(df$Vir)&!is.na(df$Res),]
df$Convergent <- df$Vir>=4 & df$Res>=2
top_sts <- names(sort(table(df$ST[df$ST!="Unknown"]),decreasing=TRUE))[1:7]
df$ST_group <- ifelse(df$ST %in% top_sts, df$ST, "Other")
bubble_df <- df %>% group_by(Res,Vir,Convergent,ST_group) %>% summarise(count=n(),.groups="drop")
wong_pal <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
st_cols <- setNames(c(wong_pal,"grey70"),c(top_sts,"Other"))
fig7 <- ggplot() +
  annotate("rect",xmin=1.7,xmax=3.3,ymin=3.7,ymax=5.3,fill="#FDEDEC",colour="#C0392B",linetype="dashed",linewidth=0.8,alpha=0.5) +
  annotate("text",x=2.5,y=5.22,label="Convergent zone [n=20]",colour="#C0392B",fontface="bold",size=6,vjust=0) +
  geom_point(data=filter(bubble_df,!Convergent),aes(x=Res,y=Vir,size=count,fill=ST_group),shape=21,colour="grey40",alpha=0.85,stroke=0.7) +
  geom_point(data=filter(bubble_df,Convergent),aes(x=Res,y=Vir,size=count,fill=ST_group),shape=24,colour="#922B21",alpha=0.95,stroke=1.2) +
  scale_fill_manual(values=st_cols,name="Sequence Type",guide=guide_legend(override.aes=list(shape=21,size=5))) +
  scale_size_area(max_size=22,name="No. isolates",breaks=c(1,5,10,20,40),guide=guide_legend(override.aes=list(shape=21,colour="grey40"))) +
  scale_x_continuous(breaks=0:3,limits=c(-0.5,3.6),labels=c("0","1","2","3")) +
  scale_y_continuous(breaks=0:5,limits=c(-0.7,5.5)) +
  annotate("text",x=-0.45,y=-0.65,label="Circle = Non-convergent     Triangle = Convergent",size=5,colour="grey40",hjust=0) +
  labs(x="Resistance Score (0-3)",y="Virulence Score (0-5)") +
  theme_classic(base_size=17) +
  theme(axis.title=element_text(face="bold",size=16),axis.text=element_text(size=15,colour="grey20"),panel.grid.major=element_line(colour="grey92",linewidth=0.4),legend.position="right",legend.title=element_text(face="bold",size=14),legend.text=element_text(size=13),plot.margin=margin(10,15,10,10))
ggsave(file.path(OUT,"Fig7_convergence.png"),fig7,width=12,height=9,dpi=300,bg="white",type="cairo")
message("Fig7 saved")
