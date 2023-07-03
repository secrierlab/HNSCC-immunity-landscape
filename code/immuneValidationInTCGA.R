library(TCGA2STAT)
library("RColorBrewer")
library(ComplexHeatmap)
library(GSVA)
library(gdata)
library(pheatmap)
library(ConsensusTME)
library(ggpubr)
library("survminer")
library("survival")
library(TCGAbiolinks)

rnaseq.hn <- getTCGA(disease="HNSC", data.type="RNASeq2", 
                     type="TPM", clinical = TRUE)
expr.hn <- t(log2(rnaseq.hn$dat+1))
clin.hn <- rnaseq.hn$clinical
clin.hn <- data.frame(clin.hn)
clin.hn$Sample <- rownames(clin.hn)

expr.hn <- data.frame(expr.hn)
dat.tcga.expr <- expr.hn[which(sapply(rownames(expr.hn), 
                                   function(x) ifelse(substr(strsplit(as.character(x),"-")[[1]][4],1,2)=="01",
                                                      "keep","discard"))=="keep"),]

query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)
data <- GDCprepare(query = query)
expr_hn <- data.frame(data@assays@data)
expr <- t(expr_hn[,-c(1:2)])
colnames(expr) <- data@rowRanges$gene_id
rownames(expr) <- data@colData$patient
expr <- data.frame(expr)
expr$Sample <- rownames(expr)

ctme <- ConsensusTME::consensusTMEAnalysis(t(as.matrix(dat.tcga.expr)), 
                                           cancer = "HNSC", 
                                           statMethod = "ssgsea")
toremove <- c("Eosinophils","Mast_cells","Neutrophils",
              "Endothelial","Fibroblasts")

# keep only categories that have sufficient coverage:
ctme <- ctme[which(!(rownames(ctme) %in% toremove)),]

save(ctme, file="ctme.TCGA.log2.RData")

pdf("plots2020/validationTCGA.heatmap.log2.TPM.pdf",w=12,h=9)
p <- pheatmap(ctme, show_colnames = FALSE)
dev.off()

ctme <- ctme[which(!(rownames(ctme) %in% c("Macrophages","Monocytes"))),]

library(ConsensusClusterPlus)
ctme.robust = ConsensusClusterPlus(as.matrix(ctme),maxK=5,reps=1000,pItem=0.80,pFeature=0.80,
                                   title="CTME robust TCGA",clusterAlg="hc",distance="euclidean",
                                   finalLinkage="complete",innerLinkage="complete",
                                   seed=1262118388.71279,plot="png")

save(ctme.robust, file="ctme.robust.log2.RData")

hr <- ctme.robust[[5]][["consensusTree"]] # of 5 because I want 5 clusters
hc <- hclust(as.dist(1-cor(t(ctme.robust), method="pearson")), method="complete")  
mycl <- cutree(hr, k=5);# h=max(hr$height)/1.5); 

# the groups are assigned according to the clusters defined by consensus clustering:
dat.annot <- data.frame(Sample = sapply(colnames(ctme),
                                        function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-")), 
                        Group = mycl)
dat.annot$SampleID <- colnames(ctme)
expr$sample <- sapply(expr$Sample, function(x) paste(strsplit(x,"\\.")[[1]][1:3],collapse="-"))
dat.annot <- merge(unique(dat.annot[,c("Sample","Group")]), expr[,c("sample","CD274")],
                   by.x="Sample",by.y="sample",
                   all.x=FALSE,all.y=FALSE)
dat.annot <- data.frame(dat.annot)
dat.annot <- dat.annot[!duplicated(dat.annot$Sample),]
rownames(dat.annot) <- dat.annot$Sample
dat.annot <- dat.annot[,-c(1)]
dat.annot$CD274 <- log2(dat.annot$CD274+1)
save(dat.annot, file="dat.annot.TCGA.log2.RData")

load("dat.annot.TCGA.log2.RData")
dat.tcga.expr <- data.frame(dat.tcga.expr)
dat.tcga.expr$Sample <- sapply(rownames(dat.tcga.expr), function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
dat.annot <- data.frame(dat.annot)
dat.annot$Sample <- rownames(dat.annot)


# library(NMF)
# ##### Plot CTME clusters robustly:
# pdf("plots2020/TCGA.ctme.robust.RSEM.pdf", onefile = FALSE)
# annClassif <- dat.annot
# hm <- aheatmap(as.matrix(ctme),
#                Rowv=as.dendrogram(hc), 
#                Colv=as.dendrogram(hr),
#                scale="none",
#                labCol = NA,
#                fontsize=5, cellheight=15,
#                cexCol = 0.9, annCol = annClassif)
# dev.off()


library(ComplexHeatmap)
library(circlize)
colnames(dat.annot) <- c("Group","PDL1")
dat.annot <- data.frame(dat.annot)
colours <- list("Group"=c("1"="#4A6FA5", 
                                  "2"="#DCD6F7",
                                  "3"="#F98948",
                                  "4"="#993955",
                          "5"="yellow"),
                "PDL1"=colorRamp2(c(0, max(dat.annot$PDL1,na.rm=T)),
                                   c("#EEE2DF", "#55251D"))
                )
ha <- HeatmapAnnotation(
  df=dat.annot,col=colours,
  annotation_name_side = "left")
ht <- Heatmap(as.matrix(ctme), name = "Infiltration",
              cluster_rows = as.dendrogram(hc),
              cluster_columns=as.dendrogram(hr),
              show_column_names=FALSE,
              top_annotation =ha,
              na_col = "yellow")
pdf("plots2021.extra/TCGA.heatmap.annot.log2.pdf",w=10,h=6)
draw(ht)
dev.off()

library(ggplot2)
library(ggpubr)
dat.annot$ImmunityGroup <- sapply(dat.annot$Group,
                                  function(x) ifelse(x==1,"MacrophageRich",
                                                     ifelse(x==2,"Absent",
                                                            ifelse(x==3,"DCrich",ifelse(x==4,"PCrich","Cytotoxic")))))
dat.annot$ImmunityGroup <- factor(dat.annot$ImmunityGroup,
                                  levels=c("Absent","MacrophageRich","DCrich","PCrich","Cytotoxic"))

my_comparisons4 <- list(c("Absent","MacrophageRich"),c("MacrophageRich","DCrich"),
                        c("DCrich","PCrich"),c("PCrich","Cytotoxic"))
pdf("plots2021.extra/TCGA.comparingPDL1status.pdf",w=6,h=6)
ggboxplot(dat.annot, 
          x = "ImmunityGroup", 
          y = "PDL1",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill="ImmunityGroup")+
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+
  xlab("")+
  ylab("PD-L1 expression (log2 TPM)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y=0,label.x= 1)
dev.off()

dat.annot.egfr <- merge(dat.annot, dat.tcga.expr[,c("Sample","EGFR")],
                        by.x="Sample",by.y="Sample",
                        all.x=FALSE,all.y=FALSE)

pdf("plots.revision/TCGA.comparingEGFRexpression.pdf",w=6,h=6)
ggboxplot(dat.annot.egfr, 
          x = "ImmunityGroup", 
          y = "EGFR",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill="ImmunityGroup")+
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+
  xlab("")+
  ylab("EGFR expression (log2 TPM)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y=0,label.x= 1)
dev.off()


## Finally, assess TLS signatures in TCGA:

# load TCGA expression data:
load("../scripts_forpub/TCGA.dat.tcga.expr.RData")

# From Helen
TLS.tfh <- c("CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1")
TLSb.th1 <- c("CD4","CCR5","CXCR3","CSF2","IGSF6","IL2RA","CD38","CD40","CD5",
              "MS4A1","SDC1","GFI1","IL1R1","IL1R2","IL10","CCL20","IRF4","TRAF6","STAT5A")
TLSchemokine <- c("CXCL13","CXCL11","CXCL10","CXCL9","CCL21","CCL19","CCL18","CCL8",
                  "CCL5","CCL4","CCL3","CCL2")


tcga.tls <- dat.tcga.expr[,unique(c(TLS.tfh,TLSb.th1,TLSchemokine))]
tcga.tls <- data.frame(tcga.tls)
tcga.tls$TLS.tfh <- apply(tcga.tls[,TLS.tfh], 1,
                               function(x) mean(x))
tcga.tls$TLSb.th1 <- apply(tcga.tls[,TLS.tfh], 1,
                          function(x) mean(x))
tcga.tls$TLSchemokine <- apply(tcga.tls[,TLSchemokine], 1,
                               function(x) mean(x))
tcga.tls$Sample <- sapply(rownames(tcga.tls), function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))

dat.annot$Sample <- rownames(dat.annot)

dat.tcga.tls <- merge(dat.annot, tcga.tls[,c("Sample",
                                             "TLS.tfh",
                                             "TLSb.th1",
                                             "TLSchemokine")],
                      by.x="Sample", by.y="Sample",
                      all.x=FALSE, all.y=FALSE)


pdf("plots2021.extra/TCGA.comparingTLS.tfh.pdf",w=6,h=6)
ggboxplot(dat.tcga.tls, 
          x = "ImmunityGroup", 
          y = "TLS.tfh",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill="ImmunityGroup")+
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+
  xlab("")+
  ylab("TLS TFH signature score")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y=0,label.x= 1)
dev.off()

pdf("plots2021.extra/TCGA.comparingTLSb.th1.pdf",w=6,h=6)
ggboxplot(dat.tcga.tls, 
          x = "ImmunityGroup", 
          y = "TLSb.th1",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill="ImmunityGroup")+
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+
  xlab("")+
  ylab("TLS TH1 signature score")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y=0,label.x= 1)
dev.off()

pdf("plots2021.extra/TCGA.comparingTLSchemokine.pdf",w=6,h=6)
ggboxplot(dat.tcga.tls, 
          x = "ImmunityGroup", 
          y = "TLSchemokine",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill="ImmunityGroup")+
  stat_compare_means(comparisons = my_comparisons4, label = "p.signif")+
  xlab("")+
  ylab("TLS chemokine score")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  stat_compare_means(label.y=0,label.x= 1)
dev.off()

### Survival:
load("../scripts_forpub/TCGA.clin.hn.RData")
clin.hn$OS <- as.numeric(apply(clin.hn[,c("daystodeath","daystolastfollowup")],
                    1, function(x) ifelse(is.na(x[1]),x[2],x[1])))
clin.hn$vitalstatus <- as.numeric(as.character(clin.hn$vitalstatus))
dat.annot$Sample <- rownames(dat.annot)
dat.annot.surv <- merge(dat.annot, clin.hn[,c("Sample","OS","vitalstatus")],
                        by.x="Sample",by.y="Sample",
                        all.x=FALSE, all.y=FALSE)
dat.annot.surv$Immunity <- sapply(dat.annot$ImmunityGroup,
                                  function(x) ifelse(x %in% c("Absent","MacrophageRich"),"Low","High"))

library(survminer)
library(survival)
fit <- survfit(Surv(OS, vitalstatus) ~ ImmunityGroup, data=dat.annot.surv)
pdf("plots2021.extra//TCGAsurvival.5groups.pdf",onefile = FALSE)
ggsurvplot(fit, linetype="strata",
           palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                       "#F98948","#993955"),
           #conf.int = TRUE, 
           pval=TRUE,#legend.labs=c("Hi", "Sex=2"),
           ggtheme = theme_minimal(), risk.table=TRUE, risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           #conf.int.style = "step",  # customize style of confidence intervals
           surv.median.line = "hv",  # add the median survival pointer.
           legend.labs = 
             c("Absent", "Macrophage rich", "DC rich", "PC rich", "Cytotoxic"))
dev.off()

fit2 <- survfit(Surv(OS, vitalstatus) ~ Immunity, data=dat.annot.surv)
pdf("plots2021.extra/TCGAsurvival.2groups.pdf",onefile = FALSE)
ggsurvplot(fit2, linetype="strata",
           palette = c("#993955", "#336699"),
           #conf.int = TRUE, 
           pval=TRUE,#legend.labs=c("Hi", "Sex=2"),
           ggtheme = theme_minimal(), risk.table=TRUE, risk.table.y.text.col = TRUE,
           risk.table.height = 0.25, # the height of the risk table
           risk.table.y.text = FALSE,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = TRUE,      # plot the number of censored subjects at time t
           ncensor.plot.height = 0.25,
           #conf.int.style = "step",  # customize style of confidence intervals
           surv.median.line = "hv"  # add the median survival pointer.
           )
dev.off()


