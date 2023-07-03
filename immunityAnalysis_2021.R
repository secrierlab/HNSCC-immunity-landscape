library(ggplot2)
library(reshape)
library(ggpubr)
library(forcats)
library(corrplot)

# Load data with CTME annotation:
load("dat.annot.plusimm.2021.RData")

# Extract immune infiltrate values:
tab.imm <- dat.annot.plusimm[,c("B_cells","Cytotoxic_cells",
                                "Dendritic_cells",#"Macrophages",
                                "NK_cells","T_cells_CD4",
                                "T_cells_CD8","T_cells_gamma_delta",
                                "T_regulatory_cells","Macrophages_M1",
                                "Macrophages_M2",#"Monocytes",
                                "Plasma_cells")]

## Dimensionality reduction did not yield anything interesting:
# # Try PCA:
# pc <- prcomp(tab.imm)
# pdf("plots2021/pca.pdf")
# plot(pc$x[,1],pc$x[,2])
# dev.off()
# 
# # Try UMAP:
# library(umap)
# library(M3C)
# pdf("plots2021/umap.pdf")
# umap(t(tab.imm))
# dev.off()
# 
# # Try tSNE:
# library(ggplot2)
# pdf("plots2021/tsne.pdf")
# tsne(t(tab.imm),perplex=20)
# dev.off()
# 
# library(phateR)
# dat.phate <- phate(tab.imm)

dat.annot.plusimm$ImmunityGroup <- factor(dat.annot.plusimm$ImmunityGroup,
                            levels=c("Absent","MacrophageRich","DCrich","PCrich","Cytotoxic"))

### Compare various parameters:
# dat.imm.red <- dat.annot.plusimm[,c(1:8,10)]
# df.melt <- melt(dat.imm.red)
# 
# newlabs <- c("Productive clonality",
#              "Observed richness",
#              "PD-L1 staining")
# names(newlabs) <- unique(df.melt$variable)

my_comparisons <- list(c("Absent","MacrophageRich"),
                       c("MacrophageRich","DCrich"),
                       c("DCrich","PCrich"),
                       c("PCrich","Cytotoxic"))

load("bulkExpMatrix.162.RData")
dat.hn <- data.frame(t(bulkExpMatrix))
dat.hn$Sample <- sapply(rownames(dat.hn),
                        function(x) paste(strsplit(x,"\\.")[[1]],collapse="-"))

dat.merged <- merge(dat.hn, 
                    dat.annot.plusimm,
                    all.x=FALSE, all.y=FALSE,
                    by.x="Sample", by.y="Sample")
df.merged.melt1 <- melt(dat.merged[,c("Sample",
                                     "TBX21","EOMES","TCF7",
                                     "ImmunityGroup")])

pdf("plots2021/compareImmG.TBETEOMES.pdf",
    w=10,h=6)
ggboxplot(df.merged.melt1, x = "ImmunityGroup", 
          y = "value",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill = "ImmunityGroup")+
  facet_wrap(~variable,scales = "free")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Expression")+
  theme(axis.text.x=element_blank(),
        axis.ticks = element_blank())+
  stat_compare_means(label.y = -1, label.x= 2)
dev.off()

## And altogether:

df.merged.melt2 <- melt(dat.merged[,c("Sample",
                                      "Productive.Clonality",
                                      "observed_richness",
                                      "PD.L1.SP263.TC",
                                      "TBX21","EOMES","TCF7",
                                      "ImmunityGroup")])
newlabs <- c("Productive clonality",
             "Observed richness",
             "PD-L1 staining",
             "TBET","EOMES","TCF7")
names(newlabs) <- c("Productive.Clonality",
                    "observed_richness",
                    "PD.L1.SP263.TC",
                    "TBX21","EOMES","TCF7")

pdf("plots2021/compareImmG.allExhaustion.pdf",
    w=12,h=5)
ggboxplot(df.merged.melt2, x = "ImmunityGroup", 
          y = "value",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill = "ImmunityGroup")+
  facet_wrap(~variable,scales = "free",nrow=1,
             labeller = labeller(variable=newlabs))+
  xlab("")+
  ylab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  stat_compare_means(label.y = 0, label.x= 0)
dev.off()




# Dot plot of enrichment:
df.totest <- dat.annot.plusimm[,c("Stage","HPV.category",
                     "Smoking.Category","Tumor.Site",
                     "Sex","ImmunityGroup")]

pvals <- NULL
odds <- NULL
pairs <- NULL
features <- NULL
for (var in 1:5) {
  tb <- table(df.totest[,c(6,var)])
  features <- c(features,colnames(tb))
  for (i in 1:nrow(tb)) {
    for (j in 1:ncol(tb)) {
        pairs <- c(pairs,paste0(rownames(tb)[i],":",colnames(tb)[j]))
        f <- fisher.test(matrix(c(tb[i,j],sum(tb[i,-j]),
                           sum(tb[-i,j]),sum(tb[-i,-j])),ncol=2))
        pvals <- c(pvals,f$p.value)
        odds <- c(odds, f$estimate)

      }
    }
              
}
names(pvals) <- pairs
names(odds) <- pairs
pvals.adj <- p.adjust(pvals, method="BH")

pvals[which(pvals<0.1)]
pvals.adj[which(pvals.adj<0.1)]

df.dot <- data.frame(Pair = pairs,
                     P.adj = pvals,
                     OR = odds)
df.dot$Group <- sapply(df.dot$Pair,
                       function(x) strsplit(as.character(x),"\\:")[[1]][1])
df.dot$Feature <- sapply(df.dot$Pair,
                       function(x) strsplit(as.character(x),"\\:")[[1]][2])

df.dot$Feature <- factor(df.dot$Feature, 
                         levels = rev(c("I","II","III","IV","Unknown tumour stage",
                                    "Positive","Negative",
                                    "Never","Past","Present",
                                    "Mouth and tongue","Oropharynx",
                                    "Hypopharynx","Larynx",
                                    "Other sites","Male","Female"
                                    )))
colnames(df.dot)[2] <- "P.value" 
#df.dot[which(df.dot$Group == "PlasmaCellHigh"),]$Group<-"PC high"

df.dot$Group <- factor(df.dot$Group, 
                       levels = c("Absent",
                                  "MacrophageRich",
                                  "DCrich",
                                  "PCrich",
                                  "Cytotoxic"))
df.dot$`Odds ratio(log 2)` <- abs(log2(df.dot$OR))
df.dot[which(is.infinite(df.dot$`Odds ratio(log 2)`)),]$`Odds ratio(log 2)` <- 5

pdf("plots2021/dotplot.featureEnrichment.pdf",
    w=12,h=5)
ggplot(df.dot, aes(x = log2(OR), y = Feature)) + 
  geom_point(aes(size = `Odds ratio(log 2)`, color = P.value)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red",high ="grey") +
  ylab(NULL) +
  xlab("Odds ratio (log2)")+
  xlim(c(-6,6))+
  facet_grid(.~Group)+
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "black", size=0.5)+
  theme(strip.text.x = element_text(size = 12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=12)) 
dev.off()

#################################################
### Next, compare TLS levels between groups:

# TLS signature genes:
TLS1 <- c("CCL19","CCL21","CXCL13","CCR7",
          "SELL","LAMP3","CXCR4","CD86","BCL6")
TLS2 <- c("CD79B","CD1D","CCR6","LAT","SKAP1",
          "CETP","EIF1AY","RBP5","PTGDS","MS4A1")

# How many genes are in my expression dataset?
length(which(colnames(dat.merged) %in% TLS1))
# all of them :)
TLS2.short <- colnames(dat.merged)[which(colnames(dat.merged) %in% TLS2)]
# 7 out of 10

dat.merged$TLSsignature.1 <- apply(dat.merged[,TLS1],
                                   1, function(x) mean(x))
dat.merged$TLSsignature.2 <- apply(dat.merged[,TLS2.short],
                                   1, function(x) mean(x))

# Now compare TLS levels between groups:
df.melt2 <- melt(dat.merged[,c("ImmunityGroup",
                               "TLSsignature.1",
                               "TLSsignature.2",
                               "CD8A","CD8B")])

pdf("plots2021/compareImmG.TLS.pdf",
    w=6,h=8)
ggboxplot(df.melt2, x = "ImmunityGroup", 
          y = "value",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          fill = "ImmunityGroup")+
  facet_wrap(~variable,scales = "free")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("Immunity group")+
  ylab("Expression")+
  theme(axis.text.x=element_blank())+
  stat_compare_means(label.y = -1, label.x= 2)
dev.off()

#############################
## Compare immune cell infiltrates between sites:

df.melt3 <- melt(dat.annot.plusimm)

df.melt.imm <- melt(dat.annot.plusimm[,-c(8:14)])
df.melt.imm$Tumor.Site <- factor(df.melt.imm$Tumor.Site,
                                 levels=c("Mouth and tongue",
                                          "Larynx",
                                              "Oropharynx",
                                          "Hypopharynx",
                                              
                                              "Other sites"))

my_comparisonsbs <- list(c("Mouth and tongue","Larynx"),
                         c("Larynx","Oropharynx"),
                         c("Oropharynx","Hypopharynx"),
                         c("Hypopharynx","Other sites"))
pdf("plots2021/compareInfiltr.byTumourSite.pdf",w=12,h=12)
ggboxplot(df.melt.imm, 
         x = "Tumor.Site", 
         y = "value",
         #palette = c("#00AFBB", "#E7B800"),
         color = "Tumor.Site",
         add = "jitter")+
  facet_wrap(~variable,ncol=5,scale="free")+
  stat_compare_means(comparisons = my_comparisonsbs, label = "p.signif")+
  xlab("Site of origin")+
  ylab("Immune infiltration")+
  theme(axis.text.x=element_blank())+
  stat_compare_means(label.y = -0.3)
dev.off()

df.melt.imm$ImmunityGroup <- factor(df.melt.imm$ImmunityGroup, 
                                          levels = c("Absent",
                                                     "MacrophageRich",
                                                     "DCrich",
                                                     "PCrich",
                                                     "Cytotoxic"))

pdf("plots2021/compareInfiltr.byImmunityGroup.pdf",w=10,h=12)
ggboxplot(df.melt.imm, 
          x = "ImmunityGroup", 
          y = "value",
          palette = c("#4A6FA5", "#DCD6F7","#FFE89F",
                      "#F98948","#993955"),
          color = "ImmunityGroup",
          add = "jitter")+
  facet_wrap(~variable,ncol=5,scale="free")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Immune infiltration")+
  theme(axis.text.x=element_blank())
dev.off()



my_comparisons3 <- list(c("Positive","Negative"))
pdf("plots2021/compareInfiltr.byHPV.pdf",w=5,h=5)
ggviolin(df.melt.imm[which(df.melt.imm$variable %in%
                          c("Macrophages_M2",
                            "Plasma_cells","Dendritic_cells")),], 
         x = "HPV.category", 
          y = "value",
          palette = c("#9A98B5", "#A80874"),
          fill = "HPV.category",
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable)+
  stat_compare_means(comparisons = my_comparisons3, label = "p.signif")+
  xlab("HPV status")+
  ylab("")+
  theme(axis.text.x=element_blank())+
  stat_compare_means(label.y = 0.7)
dev.off()

#### Corrplot of immune cells:

res1 <- cor.mtest(tab.imm, conf.level = .95)
M <- cor(tab.imm)
pdf("plots2021/immuneCorrplot.pdf")
corrplot(M, p.mat = res1$p, method="color",
         type="upper",
         insig = "label_sig",
         sig.level = c(.001, .01, .05), 
         pch.cex = .9,
         pch.col = "white", 
         order = "AOE",
         tl.col = "black")
dev.off()
                  
# Some negative correlation between plasma cells and macrophages, but not signficant:
cor.test(tab.imm$Macrophages_M2,
         tab.imm$Plasma_cells)  


########################################
### Enrichment of mutations between groups:

library(gdata)
snvs <- read.delim("../scripts_forpub/vardict.pass_20200617.tsv")
table(snvs$Functional_Class)
hist(snvs$AlleleFreq)
quantile(snvs$AlleleFreq, na.rm=TRUE)

# Read sample info:
sample.info <- read.xls("../scchn_wes_nanostring.xlsx")
snvs <- merge(snvs, sample.info[,c("SampleID","SampleName")],
              all.x=FALSE, all.y = FALSE,
              by.x="Sample", by.y="SampleName")
save(snvs, file="snvs.RData")

# remove synonymous variants:
snvs.nonsyn <- snvs[which(snvs$Functional_Class != "synonymous_variant"),]

# read in cancer genes:
cancergenes <- read.delim("Census_allWed Jul 21 12_38_18 2021.tsv")
genes <-  cancergenes$Gene.Symbol

snvs.genes <- unique(snvs.nonsyn[which(snvs.nonsyn$Gene %in% genes),
                                 c("SampleID","Gene")])

df.genes <- merge(dat.merged, snvs.genes,
                  by.x="Sample", by.y="SampleID",
                  all.x=FALSE, all.y=FALSE)

genestocheck<- names(table(df.genes$Gene)[which(table(df.genes$Gene)>5/100*nrow(dat.merged))])

pvals <- NULL
odds <- NULL
names <- NULL
j <- 0
for (g in unique(genestocheck)) {
  j<-j+1
  print(j)
  for (i in unique(df.genes$ImmunityGroup)) {
    
    names <- c(names,paste0(g,":",i))

    m <- matrix(c(length(unique(df.genes[which((df.genes$Gene == g) &
                           (df.genes$ImmunityGroup == i)),]$Sample.ID)),
    length(unique(df.genes[which((df.genes$Gene == g) &
                    (df.genes$ImmunityGroup != i)),]$Sample.ID)),
    length(unique(df.genes[which((df.genes$Gene != g) &
                (df.genes$ImmunityGroup == i)),]$Sample.ID)),
    length(unique(df.genes[which((df.genes$Gene != g) &
                          (df.genes$ImmunityGroup != i)),]$Sample.ID))),
    ncol=2)
    f <- fisher.test(m)
    pvals <- c(pvals,f$p.value)
    odds <- c(pvals,f$estimate)
  }
}
names(pvals) <- names
names(odds) <- names
pvals.adj <- p.adjust(pvals, method="BH")
pvals.adj[which(pvals.adj<0.1)]
# no enrichment of drivers in any of the groups

### Next, try by hot vs cold

df.genes$ImmunityGroupHotCold <- sapply(df.genes$ImmunityGroup,
                                        function(x) 
                                          ifelse(x %in% c("Absent","DCrich","MacrophageRich"),"Cold","Hot"))

pvals <- NULL
odds <- NULL
names <- NULL
j <- 0
for (g in unique(genestocheck)) {
  j<-j+1
  print(j)
  for (i in unique(df.genes$ImmunityGroupHotCold)) {
    
    names <- c(names,paste0(g,":",i))
    
    m <- matrix(c(length(unique(df.genes[which((df.genes$Gene == g) &
                                                 (df.genes$ImmunityGroup == i)),]$Sample.ID)),
                  length(unique(df.genes[which((df.genes$Gene == g) &
                                                 (df.genes$ImmunityGroup != i)),]$Sample.ID)),
                  length(unique(df.genes[which((df.genes$Gene != g) &
                                                 (df.genes$ImmunityGroup == i)),]$Sample.ID)),
                  length(unique(df.genes[which((df.genes$Gene != g) &
                                                 (df.genes$ImmunityGroup != i)),]$Sample.ID))),
                ncol=2)
    f <- fisher.test(m)
    pvals <- c(pvals,f$p.value)
    odds <- c(pvals,f$estimate)
  }
}
names(pvals) <- names
names(odds) <- names
pvals.adj <- p.adjust(pvals, method="BH")
pvals.adj[which(pvals.adj<0.1)]

# no enrichment of drivers in any of the groups

