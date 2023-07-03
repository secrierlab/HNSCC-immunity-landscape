###################
### This script performs TLS analyis along with other more general immune features.

library(reshape)
library(ggplot2)
library(ggpubr)
library(gdata)

# From Helen
TLS.tfh <- c("CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1")
TLSb.th1 <- c("CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40","CD5",
              "MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20","IRF4","TRAF6","STAT5A")
TLSchemokine <- c("CXCL13","CXCL11","CXCL10","CXCL9","CCL21","CCL19","CCL18","CCL8",
                  "CCL5","CCL4","CCL3","CCL2")


# Read in data:
load("dat.annot.plusimm.2021.RData")
dat.annot.plusimm$ImmunityGroup <- factor(dat.annot.plusimm$ImmunityGroup,
                                          levels=c("Absent","MacrophageRich","DCrich","PCrich","Cytotoxic"))

load("bulkExpMatrix.162.RData")
dat.hn <- data.frame(t(bulkExpMatrix))
dat.hn$Sample <- sapply(rownames(dat.hn),
                        function(x) paste(strsplit(x,"\\.")[[1]],collapse="-"))
dat.merged <- merge(dat.hn, 
                    dat.annot.plusimm,
                    all.x=FALSE, all.y=FALSE,
                    by.x="Sample", by.y="Sample")

# all generic TLS genes a and b included
setdiff(TLSchemokine, colnames(dat.merged))
# a few genes missing, but not enough to affect the signature

#Additional genes, from https://www.nature.com/articles/s41586-019-1914-8:
auxiliary <- c("CD20","CD8A","CD8B","CD4","CD21","FOXP3","GITR")
aux <- auxiliary[which(auxiliary %in% colnames(dat.merged))]

TLSchemokine[which(TLSchemokine %in% colnames(dat.merged))]
# all of them

TLSa <- TLS.tfh[which(TLS.tfh %in% colnames(dat.merged))]
TLSb <- TLSb.th1[which(TLSb.th1 %in% colnames(dat.merged))]
dat.merged$TLSsignatureA <- apply(dat.merged[,TLSa],
                                   1, function(x) mean(x))
dat.merged$TLSsignatureB <- apply(dat.merged[,TLSb],
                                  1, function(x) mean(x))
dat.merged$TLSchemokine <- apply(dat.merged[,TLSchemokine],
                                  1, function(x) mean(x))
dat.merged$TLSsignature <- apply(dat.merged[,unique(c(TLSa,TLSb))],
                                 1, function(x) mean(x))

## Other IHC markers:
ihc <- read.xls("../ts_io_head_neck_clinical_data.xlsx")
dat.merged <- merge(dat.merged, ihc[,c("Sample.ID","CD8.Ct.Cells.Mm2","CD8.Im.Cells.Mm2",
                                       "CD3.Ct.Cells.Mm2","CD3.Im.Cells.Mm2",
                                       "Density.FOXP3.TC","Density.FOXP3.IM",
                                       "Density.GITR.TC","Density.GITR.IM",
                                       "ATM.Positive.10.Percent","ATM.Total.Percentage",
                                         "CCNE1.Nuclear.H.Score","CCNE1.Total.Percentage",
                                       "c.myc.cytoplasmic.H.Score","C.Myc.Total.Cytoplasmic",
                                         "LKB1.cytoplasmic.H.Score","LKB1.Total.Cytoplasmic",
                                       "PTEN.cytoplasmic.H.Score","PTEN.Cytoplasmic.Total",
                                       "Age","Age.Category","TMB","Packs.per.Year","Packs.per.Day",
                                       "Overall.Survival.months.","Censoring.Status",
                                       "PFS","Censoring.Status..PFS.")],
                    by.x="Sample", by.y="Sample.ID",
                    all.x=TRUE, all.y=FALSE)


df.melt <- melt(dat.merged[,c("ImmunityGroup",
                               "TLSsignatureA",
                               "TLSsignatureB",
                              "TLSchemokine")])

df.melt2 <- melt(dat.merged[,c("ImmunityGroup",
                              "TLSsignature",
                              "TLSchemokine")])

df.melt3 <- melt(dat.merged[,c("ImmunityGroup","CD8.Ct.Cells.Mm2","CD8.Im.Cells.Mm2",
                               "CD3.Ct.Cells.Mm2","CD3.Im.Cells.Mm2",
                               "Density.FOXP3.TC","Density.FOXP3.IM",
                               "Density.GITR.TC","Density.GITR.IM")])

df.melt5 <- melt(dat.merged[,c("ImmunityGroup",
                               aux)])

dat.merged$CCNE1.Nuclear.H.Score<- as.numeric(dat.merged$CCNE1.Nuclear.H.Score)
dat.merged$LKB1.cytoplasmic.H.Score <- as.numeric(dat.merged$LKB1.cytoplasmic.H.Score)
dat.merged$PTEN.cytoplasmic.H.Score <- as.numeric(dat.merged$PTEN.cytoplasmic.H.Score)
dat.merged$c.myc.cytoplasmic.H.Score <- as.numeric(dat.merged$c.myc.cytoplasmic.H.Score)

dat.merged$ATM.Total.Percentage <- as.numeric(dat.merged$ATM.Total.Percentage)
dat.merged$CCNE1.Total.Percentage <- as.numeric(dat.merged$CCNE1.Total.Percentage)
dat.merged$LKB1.Total.Cytoplasmic <- as.numeric(dat.merged$LKB1.Total.Cytoplasmic)
dat.merged$PTEN.Cytoplasmic.Total <- as.numeric(dat.merged$PTEN.Cytoplasmic.Total)
dat.merged$C.Myc.Total.Cytoplasmic <- as.numeric(dat.merged$C.Myc.Total.Cytoplasmic)

df.melt6 <- melt(dat.merged[,c("ImmunityGroup",
                               "CCNE1.Nuclear.H.Score",
                               "LKB1.cytoplasmic.H.Score",
                               "PTEN.cytoplasmic.H.Score",
                               "c.myc.cytoplasmic.H.Score")],
                 id.vars = "ImmunityGroup")

df.melt7 <- melt(dat.merged[,c("ImmunityGroup",
                               "ATM.Total.Percentage",
                               "CCNE1.Total.Percentage",
                               "LKB1.Total.Cytoplasmic",
                               "PTEN.Cytoplasmic.Total",
                               "C.Myc.Total.Cytoplasmic")],
                 id.vars = "ImmunityGroup")

my_comparisons <- list(c("Absent","MacrophageRich"),
                       c("MacrophageRich","DCrich"),
                       c("DCrich","PCrich"),
                       c("PCrich","Cytotoxic"))

pdf("plots2021.tlschemo/TLSsignatures_plusAux.pdf",w=8,h=5)
ggviolin(df.melt2, x = "ImmunityGroup", 
         y = "value",
         palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                     "#F98948","#993955"),
         fill = "ImmunityGroup",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("Expression")+
  theme(axis.text.x=element_blank())+
  stat_compare_means(label.y = 0, label.x= 2)
dev.off()

pdf("plots2021.tlschemo/TLSsignatures_allsigs.pdf",w=8,h=5)
ggviolin(df.melt, x = "ImmunityGroup", 
         y = "value",
         palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                     "#F98948","#993955"),
         fill = "ImmunityGroup",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("Expression")+
  theme(axis.text.x=element_blank())+
  stat_compare_means(label.y = -1, label.x= 2)
dev.off()

df.melt3$logval <- sapply(df.melt3$value,
                          function(x) ifelse(is.na(x),NA,log10(x+1)))
pdf("plots2021.tlschemo/TLSfigure_IHCvariables.boxplot.pdf",w=14,h=4)
ggboxplot(df.melt3, 
          x = "ImmunityGroup", 
          y = "logval",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill = "ImmunityGroup")+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("IHC staining value")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )+
  stat_compare_means(label.y = 0, label.x= 2)
dev.off()

pdf("plots2021.tlschemo/CD83FOXP3expr_Signatures.pdf",w=8,h=5)
ggviolin(df.melt5, x = "ImmunityGroup", 
         y = "value",
         palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                     "#F98948","#993955"),
         fill = "ImmunityGroup",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("Expression")+
  theme(axis.text.x=element_blank())+
  stat_compare_means(label.y = -1, label.x= 2)
dev.off()


pdf("plots2021.tlschemo/IHC_Hscores.boxplot.pdf",w=10,h=8)
ggboxplot(df.melt6, 
          x = "ImmunityGroup", 
          y = "value",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill = "ImmunityGroup")+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("IHC staining value")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )+
  stat_compare_means(label.y = -1, label.x= 2)
dev.off()

pdf("plots2021.tlschemo/IHC_totalStain.boxplot.pdf",w=12,h=6)
ggboxplot(df.melt7, 
          x = "ImmunityGroup", 
          y = "value",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill = "ImmunityGroup")+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("IHC staining value")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )+
  stat_compare_means(label.y = -1, label.x= 2)
dev.off()


### Next, relate this to HPV status:

df.meltHPV <- melt(dat.merged[,c("HPV.category",
                                 "Plasma_cells",
                                 "Dendritic_cells",
                                 "Macrophages_M2")])

my_comparisons2 <- list(c("Positive","Negative"))
pdf("plots2021.tlschemo/comparingHPVstatus.pdf",w=6,h=5)
ggviolin(df.meltHPV, x = "HPV.category", 
         y = "value",
         palette = c("Positive"="#A80874",
                     "Negative"="#9A98B5"),
         fill = "HPV.category",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+
  xlab("")+
  ylab("Expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y = -0.4, label.x= 1)
dev.off()


dat.merged$PDL1.status <- sapply(dat.merged$PD.L1.SP263.TC,
                                 function(x) ifelse(is.na(x),x,ifelse(x==0,"Negative","Positive")))
df.meltPDL1 <- melt(dat.merged[,c("PDL1.status",
                                  "Productive.Clonality",
                                  "observed_richness",
                                  #"Age",
                                  "TLSsignature","TLSchemokine")])
df.meltPDL1.2 <- melt(dat.merged[,c("PDL1.status",
                                  "Productive.Clonality",
                                  "observed_richness",
                                  #"Age",
                                  "TLSsignatureA","TLSsignatureA","TLSchemokine")])

my_comparisons3 <- list(c("Positive","Negative"))
pdf("plots2021.tlschemo/comparingPDL1status.pdf",w=8,h=5)
ggviolin(df.meltPDL1[which(df.meltPDL1$PDL1.status != "Not Evaluable"),], 
         x = "PDL1.status", 
         y = "value",
         palette = c("Positive"="#55251D",
                     "Negative"="#EEE2DF",
                     "Not Evaluable"="grey"),
         fill = "PDL1.status",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif")+
  xlab("")+
  ylab("Expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.x= 1)
dev.off()

pdf("plots2021.tlschemo/comparingPDL1status.detailed.pdf",w=8,h=5)
ggviolin(df.meltPDL1.2[which(df.meltPDL1.2$PDL1.status != "Not Evaluable"),], 
         x = "PDL1.status", 
         y = "value",
         palette = c("Positive"="#55251D",
                     "Negative"="#EEE2DF",
                     "Not Evaluable"="grey"),
         fill = "PDL1.status",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif")+
  xlab("")+
  ylab("Expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.x= 1)
dev.off()


df.meltSmokers <- melt(dat.merged[,c("Smoking.Category",
                                  "observed_richness","TLSsignatureA")])

my_comparisons4 <- list(c("Never","Past"),c("Past","Present"),
                        c("Never","Present"))
pdf("plots2021.tlschemo/comparingSmokerstatus.pdf",w=6,h=6)
ggboxplot(df.meltSmokers,#[which(!is.na(df.meltSmokers$Smoker)),], 
         x = "Smoking.Category", 
         y = "value",
         palette = c("Never"="#deebf7",
                     "Past" = "#9ecae1",
                     "Present" = "#3182bd"),
         fill = "Smoking.Category")+
  facet_wrap(~variable,scales = "free",nrow=1)+
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+
  xlab("")+
  ylab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y=0,label.x= 1)
dev.off()

pdf("plots2021.tlschemo/TMBvsprodclon.pdf")
ggscatter(dat.merged[which(dat.merged$TMB<250),], 
          x = "TMB", y = "Productive.Clonality",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/TMBvsobsrichness.pdf")
ggscatter(dat.merged[which(dat.merged$TMB<250),],
          x = "TMB", y = "observed_richness",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/TMBvsobsrichness.log10.pdf")
dat.merged$observed_richness_log <- log10(dat.merged$observed_richness+1)
ggscatter(dat.merged[which(dat.merged$TMB<250),],
          x = "TMB", y = "observed_richness_log",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/Packyearsvsobsrichness.pdf")
ggscatter(dat.merged, 
          x = "Packs.per.Year", y = "observed_richness",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/PackyearsvsPRODclon.pdf")
ggscatter(dat.merged, 
          x = "Packs.per.Year", y = "Productive.Clonality",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/Packsperdaysvsobsrichness.pdf")
ggscatter(dat.merged, 
          x = "Packs.per.Day", y = "observed_richness",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/PacksperdayvsPRODclon.pdf")
ggscatter(dat.merged, 
          x = "Packs.per.Day", y = "Productive.Clonality",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/TLSvsobsrichness.pdf")
ggscatter(dat.merged[which(!is.na(dat.merged$observed_richness)),], 
          x = "TLSsignature", y = "observed_richness",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/agevsobsrichness.pdf")
ggscatter(dat.merged[which(!is.na(dat.merged$Age)),], 
          x = "Age", y = "observed_richness",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/TLSvsprodclon.pdf")
ggscatter(dat.merged[which(!is.na(dat.merged$Productive.Clonality)),], 
          x = "TLSsignature", y = "Productive.Clonality",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE) 
dev.off()

pdf("plots2021.tlschemo/TLSvsprodclon_byImmg.withlines.pdf")
ggscatter(dat.merged[which(!is.na(dat.merged$Productive.Clonality)),], 
          x = "TLSsignature", y = "Productive.Clonality",
          #add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE,
          color = "ImmunityGroup", palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                                               "#F98948","#993955"),
          label = "ImmunityGroup", repel = TRUE)+
  ylab("Productive clonality")+
  xlab("TLS signature score")+
  geom_hline(yintercept=0.085, linetype="dashed", color = "grey")+
  geom_vline(xintercept=6.91, linetype="dashed", color = "grey")
# median is 6.91
dev.off()

dat.merged1 <- dat.merged[which(!is.na(dat.merged$Productive.Clonality)),]

dat.merged1$TLSPCgroup <- NA
dat.merged1[which((dat.merged1$TLSsignature>6.91)&
                   (dat.merged1$Productive.Clonality>0.085)),]$TLSPCgroup <-"Cytotoxic"
dat.merged1[which((dat.merged1$TLSsignature<=6.91)&
                   (dat.merged1$Productive.Clonality<=0.085)),]$TLSPCgroup <-"Absent"
dat.merged1[which((dat.merged1$TLSsignature<=6.91)&
                   (dat.merged1$Productive.Clonality>0.085)),]$TLSPCgroup <-"MacrophageRich"
dat.merged1[which((dat.merged1$TLSsignature>6.91)&
                   (dat.merged1$Productive.Clonality<=0.085)),]$TLSPCgroup <-"DCrich"

dat.merged2 <- dat.merged[which(!is.na(dat.merged$observed_richness)),]


pdf("plots2021.tlschemo/TLSvsobsrich_byImmg.pdf")
ggscatter(dat.merged2, 
          x = "TLSsignature", y = "observed_richness",
          #add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE,
          color = "ImmunityGroup", palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                                               "#F98948","#993955"),
          label = "ImmunityGroup", repel = TRUE)+
  ylab("Observed richness")+
  xlab("TLS signature score")
dev.off()

dat.merged2$observed_richness_log <- log10(dat.merged2$observed_richness+1)
pdf("plots2021.tlschemo/TLSvsobsrich_byImmg.log10.pdf")
ggscatter(dat.merged2, 
          x = "TLSsignature", y = "observed_richness_log",
          #add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE,
          color = "ImmunityGroup", palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                                               "#F98948","#993955"),
          label = "ImmunityGroup", repel = TRUE)+
  ylab("Observed richness (log 10)")+
  xlab("TLS signature score")+
  geom_hline(yintercept=2.8 , linetype="dashed", color = "grey")+
  geom_vline(xintercept=6.9, linetype="dashed", color = "grey")
dev.off()

pdf("plots2021.tlschemo/TLSvsobsrich_byImmg.log10.plushist.pdf",w=10,h=10)
sp <- ggscatter(dat.merged2, 
          x = "TLSsignature", y = "observed_richness_log",
          #add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE,
          color = "ImmunityGroup", palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                                               "#F98948","#993955"),
          label = "ImmunityGroup", repel = TRUE)+
  ylab("Observed richness (log 10)")+
  xlab("TLS signature score")+
  geom_hline(yintercept=3.1, linetype="dashed", color = "grey")+
  geom_vline(xintercept=6.97, linetype="dashed", color = "grey")
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(dat.merged2, "TLSsignature", fill = "ImmunityGroup",
                   palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                               "#F98948","#993955"))
yplot <- ggdensity(dat.merged2, "observed_richness_log", fill = "ImmunityGroup", 
                   palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                               "#F98948","#993955"))+
  rotate()

# Cleaning the plots
sp <- sp
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
library(cowplot)
plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))

dev.off()


dat.merged2$TLSORgroup <- NA
dat.merged2[which((dat.merged2$TLSsignature>6.97)&
             (dat.merged2$observed_richness_log>3.1)),]$TLSORgroup <-"Cytotoxic"
dat.merged2[which((dat.merged2$TLSsignature<=6.97)&
                   (dat.merged2$observed_richness_log<=3.1)),]$TLSORgroup <-"Absent"
dat.merged2[which((dat.merged2$TLSsignature<=6.97)&
                   (dat.merged2$observed_richness_log>3.1)),]$TLSORgroup <-"Mixed"
dat.merged2[which((dat.merged2$TLSsignature>6.97)&
                   (dat.merged2$observed_richness_log<=3.1)),]$TLSORgroup <-"DCrich"

pdf("plots2021.tlschemo/TLSvsprodclon_byPDL1.pdf")
ggscatter(dat.merged[which(!is.na(dat.merged$Productive.Clonality)),], 
          x = "TLSsignature", y = "Productive.Clonality",
          #add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE,
          color = "PDL1.status", palette = c("Positive"="#55251D",
          "Negative"="#C9B7AD",
          "Not Evaluable"="grey"),
          label = "PDL1.status", repel = TRUE)+
  ylab("Productive clonality")+
  xlab("TLS signature score")
dev.off()


#Survival
require("survival")
library(survminer)

q<-quantile(dat.merged1$Productive.Clonality,na.rm=T)
dat.merged1$ProdClonhighlow <- sapply(dat.merged1$Productive.Clonality,
                                     function(x) 
                                       ifelse(x<q[3],
                                              "ClonalityLow","ClonalityHigh"))
dat.merged1$GroupComplex <- apply(dat.merged1[,c("ImmunityGroup","ProdClonhighlow")],
                                  1, function(x) ifelse(x[1] %in% c("Absent","MacrophageRich"),
                                                        paste0("ImmuneLow:",x[2]),
                                                        paste0("ImmuneHigh:",x[2])))

dat.merged.keep <- dat.merged1[which(!is.na(dat.merged1$Productive.Clonality)),]

fit<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ GroupComplex, data = dat.merged.keep)

q<-quantile(dat.merged$TLSsignature,na.rm=T)
dat.merged$TLShighlow <- sapply(dat.merged$TLSsignature,
                                     function(x) 
                                       ifelse(x<q[3],
                                              "TLSLow","TLSHigh"))
dat.merged$GroupComplexTLS <- apply(dat.merged[,c("ImmunityGroup","TLShighlow")],
                                 1, function(x) ifelse(x[1] %in% c("Absent","MacrophageRich"),
                                                       paste0("ImmuneLow:",x[2]),
                                                       paste0("ImmuneHigh:",x[2])))

fit2<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ GroupComplexTLS, data = dat.merged)

# Drawing survival curves. Change color, linetype by strata, risk.table color by strata
pdf("plots2021.tlschemo/survival.byGroupComplexTLS.pdf",onefile = FALSE,w=20)
ggsurvplot(fit2,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()

### Surv by TLS/obs rich combination:
dat.merged2$TLSORgroupbinary <- sapply(dat.merged2$TLSORgroup,
                                      function(x) ifelse(is.na(x),NA,ifelse(x=="Absent",x,"Other")))
dat.mergedsubset <- dat.merged2[which(!(dat.merged2$TLSORgroup == "Mixed")),]
dat.mergedsubset$TLSORgroupbinary <- sapply(dat.mergedsubset$TLSORgroup,
                                       function(x) ifelse(is.na(x),NA,ifelse(x=="Absent",x,"Other")))
dat.mergedsubset2 <- dat.merged2[which(!(dat.merged2$TLSORgroup %in% c("Mixed","Cytotoxic"))),]


dat.merged1$TLSPCgroupbinary <- sapply(dat.merged1$TLSPCgroup,
                                      function(x) ifelse(is.na(x),NA,ifelse(x=="Absent",x,"Other")))

dat.merged$ProdClonBinary <- sapply(dat.merged$productive_clonality,
                                      function(x) ifelse(is.na(x),NA,ifelse(x>3.1,"high","low")))
dat.merged$TLSBinary <- sapply(dat.merged$TLSsignature,
                                    function(x) ifelse(is.na(x),NA,ifelse(x>6.97,"high","low")))


fit4<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSORgroup, data = dat.merged2)
fit4subset<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSORgroup, data = dat.mergedsubset)
fit4subset2<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSORgroup, data = dat.mergedsubset2)
fit4sub <- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSORgroupbinary, data = dat.merged2)
fit4sub.subset<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSORgroupbinary, data = dat.mergedsubset)
fit5<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSPCgroup, data = dat.merged1)
fit6<- survfit(Surv(Overall.Survival.months., Censoring.Status) ~ TLSPCgroupbinary, data = dat.merged1)

dat.merged2$Censoring.Status..PFS. <- sapply(dat.merged2$Censoring.Status..PFS.,
                                             function(x) ifelse(is.na(x),NA, ifelse(x=="DiseaseFree",0,1)))
dat.merged1$Censoring.PFS <- sapply(dat.merged1$Censoring.Status..PFS.,
                                             function(x) ifelse(is.na(x),NA, ifelse(x=="DiseaseFree",0,1)))
fit4pfs<- survfit(Surv(PFS, Censoring.Status..PFS.) ~ TLSORgroup, data = dat.merged2)
fit4subpfs<- survfit(Surv(PFS, Censoring.Status..PFS.) ~ TLSORgroupbinary, data = dat.merged2)
fit4subpfs.subset<- survfit(Surv(PFS, Censoring.Status..PFS.) ~ TLSORgroupbinary, data = dat.mergedsubset)
fit5pfs<- survfit(Surv(PFS, Censoring.PFS) ~ TLSPCgroup, data = dat.merged1)
fit6pfs<- survfit(Surv(PFS, Censoring.PFS) ~ TLSPCgroupbinary, data = dat.merged1)

# Drawing survival curves. Change color, linetype by strata, risk.table color by strata
pdf("plots2021.tlschemo/survival.byTLSORgroup.pdf",onefile = FALSE,w=6,h=10)
ggsurvplot(fit4,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_minimal(), 
           risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           palette = c("#4A6FA5","#FFE89F","#993955",
                       "grey"))
dev.off()

pdf("plots2021.tlschemo/survival.byTLSORgroup.OS.pdf",onefile = FALSE,w=6,h=10)
ggsurvplot(fit4,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_minimal(), 
           risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           palette = c("#4A6FA5","#FFE89F","#993955",
                       "grey"))
dev.off()

pdf("plots2021.tlschemo/survival.byTLSORgroup.OS.subset.pdf",onefile = FALSE,w=6,h=10)
ggsurvplot(fit4subset,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_minimal(), 
           risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           palette = c("#4A6FA5","#FFE89F","#993955",
                       "grey"))
dev.off()

pdf("plots2021.tlschemo/survival.byTLSORgroup.binary.pdf",onefile = FALSE,w=6,h=10)
ggsurvplot(fit4sub,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_minimal(), 
           risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           palette = c("#4A6FA5","#993955"))
dev.off()

pdf("plots2021.tlschemo/survival.byTLSORgroup.binary.OS.pdf",onefile = FALSE,w=6,h=10)
ggsurvplot(fit4sub,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_minimal(), 
           risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           palette = c("#4A6FA5","#993955"))
dev.off()

pdf("plots2021.tlschemo/pfs.byTLSORgroupBinary.pdf",onefile = FALSE,w=6,h=6)
ggsurvplot(fit4subpfs,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(),
           palette = c("#4A6FA5","#993955"))
dev.off()

pdf("plots2021.tlschemo/pfs.byTLSORgroupBinary.subset.pdf",onefile = FALSE,w=6,h=6)
ggsurvplot(fit4subpfs.subset,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(),
           palette = c("#4A6FA5","#993955"))
dev.off()

pdf("plots2021.tlschemo/pfs.byTLSORgroupBinary.subset.OS.pdf",onefile = FALSE,w=6,h=6)
ggsurvplot(fit4sub.subset,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(),
           palette = c("#4A6FA5","#993955"))
dev.off()

pdf("plots2021.tlschemo/pfs.byTLSORgroup.pdf",onefile = FALSE,w=6)
ggsurvplot(fit4pfs,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw(),
           palette = c("#4A6FA5","#FFE89F","#993955",
                       "grey"))
dev.off()


pdf("plots2021.tlschemo/survival.byTLSPCgroup.pdf",onefile = FALSE,w=6)
ggsurvplot(fit5,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()

pdf("plots2021.tlschemo/pfs.byTLSPCgroup.pdf",onefile = FALSE,w=6)
ggsurvplot(fit5pfs,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()

pdf("plots2021.tlschemo/survival.byTLSPCgroupBinary.pdf",onefile = FALSE,w=6)
ggsurvplot(fit6,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()

pdf("plots2021.tlschemo/pfs.byTLSPCgroupBinary.pdf",onefile = FALSE,w=6)
ggsurvplot(fit6pfs,          
           pval = TRUE, #conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()
