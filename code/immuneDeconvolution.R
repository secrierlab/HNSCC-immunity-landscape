library(ConsensusTME)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(stringr)
library(ConsensusClusterPlus)
library(NMF)
library(ComplexHeatmap)
library(gdata)


## Load expression data:
dat.expr <- read.csv("../SCCHN-nanostring_data.csv", stringsAsFactors = FALSE)
bulkExpMatrix <- as.matrix(dat.expr[,-1])
bulkExpMatrix <- log2(bulkExpMatrix+1)
rownames(bulkExpMatrix) <- dat.expr$Gene

missingPanel2 <- names(which(apply(bulkExpMatrix,2, function(x) sum(is.na(x)))>0))

# Keep only samples with both panels measured:

bulkExpMatrix <- bulkExpMatrix[,which(!(colnames(bulkExpMatrix) %in% missingPanel2))]
# 162 samples left
save(bulkExpMatrix, file="bulkExpMatrix.162.RData")

### Checking for overlap with gene sets in ConsensusTME:
load("../scripts_forpub/ConsensusTME-master/data/consensusGeneSets.rda")
consensusHNSC <- consensusGeneSets$HNSC
consensusHNSC$T_cells_CD8 <- c("CD8A",consensusHNSC$T_cells_CD8)

for (ct in names(consensusHNSC)) {
  print(paste("====",ct))
  int <- length(intersect(consensusHNSC[ct][[1]],rownames(bulkExpMatrix)))
  print(paste(int,"covered out of",length(consensusHNSC[ct][[1]])))
}


# Remove if more than half the genes are missing:
toremove <- c("Endothelial","Fibroblasts")
ctme <- ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "HNSC", statMethod = "ssgsea")
# Keep only categories that have sufficient coverage:
ctme <- ctme[which(!(rownames(ctme) %in% toremove)),]
# Reformat sample names:
colnames(ctme) <- sapply(colnames(ctme), function(x) str_replace_all(x,"\\.","-"))

# Load clinical annotation data:
dat.hn <- read.xls("../ts_io_head_neck_clinical_data.xlsx",
                       header=TRUE, stringsAsFactors=FALSE)
dat.older <- read.table("../HN_MASTERCOHORT_withExhGroups11-15-2017.txt", header=TRUE,
                     sep="\t", stringsAsFactors = FALSE)
dat.hn <- merge(dat.hn, dat.older[,c("Sample.ID","observed_richness")],
                by.x="Sample.ID", by.y="Sample.ID",
                all.x=TRUE, all.y=FALSE)
save(dat.hn, file="clinical.dat.hn.RData")
dat.annot <- dat.hn[which(dat.hn$Sample.ID %in% colnames(ctme)),
                    c("Sex","Race","HPV.category",
                       "Tumor.Site","Stage",
                       "Smoking.Category","Immunoscore",
                      "PD.L1.SP263.CPS","PD.L1.SP263.TC","PD.L1.SP263.IC",
                      "Productive.Clonality",
                      "observed_richness"
                       )] 
rownames(dat.annot) <- dat.hn[which(dat.hn$Sample.ID %in% colnames(ctme)),]$Sample.ID

#What if I select the top most variable cell populations?
ctme.orig <- ctme
# pops.sorted <- apply(ctme, 1, function(x) var(x))
# keep <- names(pops.sorted[pops.sorted>0.005])
# ctme <- ctme[keep,]
ctme <- ctme.orig[which(!(rownames(ctme.orig) %in% c("Macrophages"))),]

#Least variable categories:
sort(apply(ctme, 1, function(x) var(x)))
# Remove least variable categories:
ctme <- ctme[which(!(rownames(ctme) %in% c("Monocytes",
                                         "Eosinophils"))),]

# Robust clusters:
ctme.robust = ConsensusClusterPlus(as.matrix(ctme),maxK=10,reps=10000,pItem=0.80,pFeature=1,
                               title="CTME_robust",clusterAlg="hc",distance="euclidean",
                               finalLinkage="complete",innerLinkage="complete",
                               seed=1262118388.71279,plot="png")

hr <- ctme.robust[[5]][["consensusTree"]] # of 3 because I want 3 clusters
hc <- hclust(as.dist(1-cor(t(ctme), method="pearson")), method="complete")  
mycl <- cutree(hr, k=5);# h=max(hr$height)/1.5); 

# the groups are assigned according to the clusters defined by consensus clustering:
dat.annot$Group <- mycl

# Improved annot:
#colnames(dat.annot)[ncol(dat.annot)] <- "ImmunityGroup"
dat.annot$HPV.category <- factor(dat.annot$HPV.category,
                                 levels=rev(c("Positive","Negative"))) 
#dat.annot$ImmunityGroup <- factor(dat.annot$ImmunityGroup) 
dat.annot$Tumor.Site <- factor(dat.annot$Tumor.Site) 
dat.annot$Stage <- factor(dat.annot$Stage,
                                 levels=rev(c("I","II","III","IV","Unknown tumour stage"))) 
dat.annot$Immunoscore <- factor(dat.annot$Immunoscore,
                          levels=rev(c("I 0","I 1","I 2","I 3","I 4","I NaN"))) 

colours <- list(#"ImmunityGroup"=c("Absent"="#4A6FA5", 
                #                  "Weak"="#DCD6F7",
                #                  "Exhausted"="#F98948",
                 #                 "PlasmaCellHigh"="#993955"),
                "Sex" = c("Male"="#998ec3","Female"="#f1a340"),
                "Race" = c("African American"="#8dd3c7",
                           "Asian"="#bebada",
                           "Caucasian"="#ffffb3",
                           "Hispanic"="#fb8072",
                           "Unknown"="white"),
                "HPV.category"=c("Positive"="#A80874",
                              "Negative"="#9A98B5"),
                "Stage"=c("I"="#f6e8c3","II"="#d8b365",
                          "III"="#c7eae5","IV"="#5ab4ac",
                          "Unknown tumour stage"="grey"),
                "Smoking.Category"=c("Never"="#deebf7",
                           "Past" = "#9ecae1",
                           "Present" = "#3182bd"),
                "Tumor.Site"=c("Hypopharynx"="#3BB273",
                         "Larynx"="#7768AE",
                         "Mouth and tongue"="#E1BC29",
                         "Oropharynx"="#776258",
                         "Other sites"="lightgrey"),
                "Immunoscore"=c("I 0"="#67a9cf","I 1"="#d1e5f0","I 2"="#fddbc7",
                                "I 3"="#ef8a62","I 4"="#b2182b",
                                "I NaN"="white"),
                "PD.L1.SP263.CPS"=colorRamp2(c(0, max(dat.annot$PD.L1.SP263.CPS,na.rm=T)),
                                                 c("#EEE2DF", "#55251D")),
                "PD.L1.SP263.TC"=colorRamp2(c(0, max(dat.annot$PD.L1.SP263.TC,na.rm=T)),
                                             c("#EEE2DF", "#55251D")),
                "PD.L1.SP263.IC"=colorRamp2(c(0, max(dat.annot$PD.L1.SP263.IC,na.rm=T)),
                                             c("#EEE2DF", "#55251D")),
                "Productive.Clonality"=colorRamp2(c(0, max(dat.annot$Productive.Clonality,na.rm=T)),
                                                 c("#EEE2DF", "#55251D")),
                "observed_richness"=colorRamp2(c(0, max(dat.annot$observed_richness,na.rm=T)),
                                                 c("#EEE2DF", "#55251D")))

ha <- HeatmapAnnotation(
  df=dat.annot,col=colours,
  annotation_name_side = "left")
ht <- Heatmap(as.matrix(ctme), name = "Infiltration",
        cluster_rows = as.dendrogram(hc),
        cluster_columns=as.dendrogram(hr),
        show_column_names=FALSE,
        #column_split = dat.annot$HPV.category,
        top_annotation =ha,
        na_col = "yellow")
pdf("plots/consensusTME/heatmap.annotated.5clusters.2021.ssgsea.reduced.pdf",w=12,h=8)
draw(ht)
dev.off()

dat.annot$Sample <- rownames(dat.annot)

dat.hn.plusgroups <- merge(dat.hn, dat.annot[,c("Sample","Group")],
                           by.x="Sample.ID", by.y="Sample",
                           all.x=TRUE, all.y=FALSE)
dat.hn.plusgroups$Group2 <- sapply(dat.hn.plusgroups$Group,
                                   function(x) ifelse(is.na(x),x,ifelse(x %in% c(4,5),"high","low")))

dat.hn.plusgroups$ImmunityGroup <- sapply(dat.hn.plusgroups$Group,
                                           function(x) ifelse(x==1,"Absent",
                                                              ifelse(x==2,"DCrich",ifelse(x==3,"MacrophageRich",
                                                                                          ifelse(x==4,"PCrich","Cytotoxic")))))
dat.hn.plusgroups$ImmunityGroupSimple <- sapply(dat.hn.plusgroups$Group,
                                        function(x)  ifelse((x==5) | (x==2),"High","Low"))

dat.annot$ImmunityGroup <- sapply(dat.annot$Group,
                                          function(x) ifelse(x==1,"Absent",
                                                             ifelse(x==2,"DCrich",ifelse(x==3,"MacrophageRich",
                                                                                         ifelse(x==4,"PCrich","Cytotoxic")))))
dat.annot$ImmunityGroupSimple <- sapply(dat.annot$Group,
                                  function(x) ifelse((x==5) | (x==2),"High","Low"))


df.ctme <- data.frame(t(ctme))
df.ctme$Sample <- rownames(df.ctme)
dat.annot.plusimm <- merge(dat.annot,
                           df.ctme,by.x="Sample",
                           by.y="Sample",
                           all.x=FALSE, all.y=FALSE)

write.table(dat.annot.plusimm,
            file="HN.groupAnnotation.txt",sep="\t",
            quote=FALSE, row.names = FALSE)
save(dat.annot.plusimm, file="dat.annot.plusimm.2021.RData")


require("survival")
library(survminer)
dat.hn.plusgroups$AliveDead <- sapply(dat.hn.plusgroups$Survival.Status,
                                      function(x) ifelse(x=="LIVING",0,ifelse(x=="DECEASED",1,x)))
fit<- survfit(Surv(Overall.Survival.months., AliveDead) ~ ImmunityGroup, data = dat.hn.plusgroups)
fit2<- survfit(Surv(Overall.Survival.months., AliveDead) ~ Group2, data = dat.hn.plusgroups)
fit3<- survfit(Surv(Overall.Survival.months., AliveDead) ~ ImmunityGroupSimple, data = dat.hn.plusgroups)

# Drawing survival curves. Change color, linetype by strata, risk.table color by strata
pdf("plots/survival/survival.byGroupDetailed.5clusters.2021.pdf",onefile = FALSE)
ggsurvplot(fit,          
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
           #palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()

pdf("plots/survival/survival.byGroup.5clusters.2021.pdf",onefile = FALSE)
ggsurvplot(fit2,          
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
           #palette = c("#e7298a","#d95f02","#66a61e"))
dev.off()

pdf("plots/survival/survival.byGroupSimple.5clusters.2021.pdf",onefile = FALSE)
ggsurvplot(fit3,          
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#e7298a","#d95f02","#66a61e"))
dev.off()

fit4<- survfit(Surv(OS.days, AliveDead) ~ ImmunityGroup, data = dat.hn.plusgroups)
fit5<- survfit(Surv(OS.days, AliveDead) ~ ImmunityGroupSimple, data = dat.hn.plusgroups)

# Drawing survival curves. Change color, linetype by strata, risk.table color by strata
pdf("plots2020/survival.byImmunityGroup.pdf",onefile = FALSE,w=12,h=7)
ggsurvplot(fit4,          
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()

pdf("plots2020/survival.byImmunityGroupSimple.pdf",onefile = FALSE)
ggsurvplot(fit5,          
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           ggtheme = theme_bw())
#palette = c("#1b9e77","#d95f02","#66a61e","#e7298a"))
dev.off()


table(dat.annot[,c("HPV.category","ImmunityGroup")])
fisher.test(matrix(c(14,33+15,10,27+44),ncol=2))

table(dat.annot.keep[,c("Smoker","Group")])
fisher.test(matrix(c(8,15,19,32+23+29+16),ncol=2))

table(dat.annot.keep[,c("Basic.Site","Group")])
fisher.test(matrix(c(10,11,14,5+47+24+15+17),ncol=2))

table(dat.annot.keep[,c("Stage..Clinical.Simplified","Group")])
fisher.test(matrix(c(4,19,17+19,14+7+32+24),ncol=2))

# Exhaustion levels?
dat.hn.plusgroups$ImmunityGroup <- factor(dat.hn.plusgroups$ImmunityGroup,
                                          levels=c("Absent","MacrophageRich",
                                                   "PCrich","DCrich","Cytotoxic"))
my_comparisons <- list(c("Absent","MacrophageRich"),
                       c("MacrophageRich","PCrich"),
                       c("PCrich","DCrich"),
                       c("DCrich","Cytotoxic"),
                       c("DCrich","Absent"),
                       c("PCrich","Absent"))

df.bulkexpr <- data.frame(t(bulkExpMatrix))
df.bulkexpr$Sample <- sapply(rownames(df.bulkexpr),
                             function(x) str_replace_all(x,"\\.","-"))
dat.hn.plusgroups <- merge(dat.hn.plusgroups[,-c(98:911)],
                           df.bulkexpr, by.x="Sample.ID",
                           by.y="Sample",all.x=FALSE, all.y=FALSE)

### survival multivariate coxph:
model <- coxph( Surv(Overall.Survival.months., AliveDead) ~ 
                  Age+Sex+Stage+T.Category+N.Category+M.Category+
                  HPV.category+Smoking.Category+Tumor.Site+Recurrence+
                  ImmunityGroup,
                data = dat.hn.plusgroups )
pdf("plots/survival/multivariateCoxPH.pdf",onefile=FALSE,h=10)
ggforest(model)
dev.off()

pdf("plots/comp/compare.EOMES.pdf")
ggviolin(dat.hn.plusgroups, x = "ImmunityGroup", y = "EOMES", fill = "Group",
         #palette = c("#66a61e","#fb8072","#984ea3"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 9)  
dev.off()

pdf("plots/comp/compare.CD8.pdf")
ggviolin(dat.hn.plusgroups, x = "ImmunityGroup", y = "CD8.Ct.Cells.Mm2", fill = "Group",
         #palette = c("#66a61e","#fb8072","#984ea3"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 9)  
dev.off()

pdf("plots/comp/compare.CD8_IM.pdf")
ggviolin(dat.hn.plusgroups, x = "ImmunityGroup", y = "CD8.Im.Cells.Mm2", fill = "Group",
         #palette = c("#66a61e","#fb8072","#984ea3"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 9)  
dev.off()

pdf("plots2020/compare.TBET.pdf")
ggviolin(dat.hn.plusgroups, x = "ImmunityGroup", y = "TBX21", fill = "Group",
         #palette = c("#66a61e","#fb8072","#984ea3"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 9)  
dev.off()

# Load mutation data:
snvs <- read.delim("../scripts_forpub/vardict.pass_20200617.tsv")
snvs.old <- read.delim("../scripts_forpub/vardict2mut.pass.txt")
snvs.old.ids <- unique(snvs.old[,c("Sample","SampleID")])

snvs <- merge(snvs,
              snvs.old.ids,
              by.x="Sample", by.y="Sample",
              all.x=FALSE, all.y=FALSE)

snvs.tp53 <- unique(snvs[which(snvs$Gene=="TP53"),]$SampleID)
length(which(snvs.tp53 %in% dat.hn.plusgroups[which(dat.hn.plusgroups$ImmunityGroup=="Cytotoxic"),]$Sample.ID))
length(which(snvs.tp53 %in% dat.hn.plusgroups[which(dat.hn.plusgroups$ImmunityGroup=="PCrich"),]$Sample.ID))
length(which(snvs.tp53 %in% dat.hn.plusgroups[which(dat.hn.plusgroups$ImmunityGroup=="DCrich"),]$Sample.ID))
length(which(snvs.tp53 %in% dat.hn.plusgroups[which(dat.hn.plusgroups$ImmunityGroup=="Weak"),]$Sample.ID))
length(which(snvs.tp53 %in% dat.hn.plusgroups[which(dat.hn.plusgroups$ImmunityGroup=="Absent"),]$Sample.ID))

table(dat.hn.plusgroups$ImmunityGroup)

fisher.test(matrix(c(8,21-8,1+16+12,29+77+14+40-1-16-12),nrow=2))
# p-value = 0.04404
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9047687 7.9910925
# sample estimates:
#   odds ratio 
# 2.760492 
#2-fold enrichment of p53 mutations in cytotoxic group
fisher.test(matrix(c(8+1+16,21+40+14-8-1-16,
                     12,29+77-12),nrow=2))
# 4 fold enrichment in the top 3 groups
# data:  matrix(c(8 + 1 + 16, 21 + 40 + 14 - 8 - 1 - 16, 12, 29 + 77 - 12), nrow = 2)
# p-value = 0.000611
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.712323 9.258280
# sample estimates:
#   odds ratio 
# 3.885234 
fisher.test(matrix(c(8+16,21+40-8-16,
                     13,29+77+14-13),nrow=2))
# 5-fold enrichment in DCrich+cytotoxic
# data:  matrix(c(8 + 16, 21 + 40 - 8 - 16, 13, 29 + 77 + 14 - 13), nrow = 2)
# p-value = 2.278e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   2.314979 12.555781
# sample estimates:
#   odds ratio 
# 5.281002 
fisher.test(matrix(c(16,40-16,
                     13+8,21+29+77+14-13-8),nrow=2))
#4 fold enrichment in DCrich group (highest enrichment across groups)
# data:  matrix(c(16, 40 - 16, 13 + 8, 21 + 29 + 77 + 14 - 13 - 8), nrow = 2)
# p-value = 0.001383
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.597122 8.913657
# sample estimates:
#   odds ratio 
# 3.77486 

save(dat.hn.plusgroups, file="dat.hn.plusgroups.RData")
