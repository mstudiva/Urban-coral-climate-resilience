#### PACKAGES ####

# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")
# BiocManager::install("DESeq2",dependencies=T)
# BiocManager::install("arrayQualityMetrics",dependencies=T)  # requires Xquartz, xquartz.org
# BiocManager::install("BiocParallel")

# install.packages("pheatmap")
# install.packages("VennDiagram")
# install.packages("gplots")
# install.packages("vegan")
# install.packages("plotrix")
# install.packages("ape")
# install.packages("ggplot2")
# install.packages("rgl")
# install.packages("adegenet")


#### DATA IMPORT ####
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)
library(tidyverse)

#read in counts
counts = read.table("../../../../../raw/ofav/orthogroup/allcounts_host.txt")

# how many genes we have total?
nrow(counts) 
ncol(counts)

# how does the data look? 
head(counts)

keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData)
ncol(countData)
write.csv(countData, file="countData.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna)
ncol(counts4wgcna)
write.csv(counts4wgcna, file="counts4wgcna.csv")

# importing a design .csv file
design = read.csv("../../../../../raw/design_ofav.csv", head=TRUE)
design
design$site <- as.factor(design$site)
design$tank <- as.factor(design$tank)
design$ph <- as.factor(design$ph)
design$temp <- as.factor(design$temp)
design$treatment <- as.factor(design$treatment)
design$treat <- as.factor(design$treat)
str(design)


#### MODEL DESIGN and OUTLIERS ####

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+treat)

# reorders treatment factor according to "control" vs "treatment" levels
dds$site <- factor(dds$site, levels = c("Emerald", "Rainbow", "Star", "MacN"))
dds$treat <- factor(dds$treat, levels = c("CC", "LC", "CH", "LH"))

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("site"),force=T)
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested. I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# if there were outliers:
outs=c(30,46,52,59)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+treat)
dds$site <- factor(dds$site, levels = c("Emerald", "Rainbow", "Star", "MacN"))
dds$treat <- factor(dds$treat, levels = c("CC", "LC", "CH", "LH"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="initial.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,2],design[,11],sep=".")
# renames the column names
colnames(vsd)=snames

save(vsd,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ site+treat)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=TRUE)
# vsd.wg=assay(rlog(wg), blind=TRUE)
head(vsd.wg)
colnames(vsd.wg)=snames
save(vsd.wg,design,file="data4wgcna.RData")


#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_ofav_host.pdf", width=15, height=15)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
library(ape)

conditions=design
conditions$site <- factor(conditions$site, levels = c("Emerald", "Rainbow", "Star", "MacN"))
conditions$treat <- factor(conditions$treat, levels = c("CC", "LC", "CH", "LH"))

# creating a PCoA eigenvalue matrix
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors
# export this table for % variation explained by each axis (Relative_eig column)
write.csv(dds.pcoa$values, file = "PCoA_variance.csv")

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues 
pdf(file="PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's

# plotting PCoA by site and treatment
pdf(file="PCoA_ofav_host.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))],pch=c(6,2,25,17)[as.numeric(as.factor(conditions$treat))], bg= c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1 (20.0%)", ylab="Coordinate 2 (9.4%)", main="Site")
ordiellipse(scores, conditions$site, label=F, col=c("#018571","#80cdc1","#dfc27d","#a6611a"))
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), fill = c("#018571","#80cdc1","#dfc27d","#a6611a"), bty="n")
legend("topleft", legend=c("CC", "CH", "LC", "LH"), pch=c(6,2,25,17), bty="n", pt.bg = c(NA,NA,"black",NA))
plot(scores[,1], scores[,2],col=c("#92c5de","#f4a582","#0571b0","#ca0020")[as.numeric(as.factor(conditions$treat))],pch=c(15,0,1,16)[as.numeric((as.factor(conditions$site)))], xlab="Coordinate 1 (20.0%)", ylab="Coordinate 2 (9.4%)", main="Treatment")
ordiellipse(scores, conditions$treat, label=F, col=c("#92c5de","#f4a582","#0571b0","#ca0020"))
legend("topleft", legend=c("CC", "CH", "LC", "LH"), fill = c("#92c5de","#f4a582","#0571b0","#ca0020"), bty="n")
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), pch=c(15,0,1,16), bty="n")
dev.off()

# PCo's 2 and 3
pdf(file="PCoA_ofav_host_PC3.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,2], scores[,3],col=c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))],pch=c(6,2,25,17)[as.numeric(as.factor(conditions$treat))], bg= c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))], xlab="Coordinate 2 (9.4%)", ylab="Coordinate 3 (6.3%)", main="Site")
ordiellipse(scores, conditions$site, label=F, col=c("#018571","#80cdc1","#dfc27d","#a6611a"))
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), fill = c("#018571","#80cdc1","#dfc27d","#a6611a"), bty="n")
legend("topleft", legend=c("CC", "CH", "LC", "LH"), pch=c(6,2,25,17), bty="n", pt.bg = c(NA,NA,"black",NA))
plot(scores[,2], scores[,3],col=c("#92c5de","#f4a582","#0571b0","#ca0020")[as.numeric(as.factor(conditions$treat))],pch=c(15,0,1,16)[as.numeric((as.factor(conditions$site)))], xlab="Coordinate 2 (9.4%)", ylab="Coordinate 3 (6.3%)", main="Treatment")
ordiellipse(scores, conditions$treat, label=F, col=c("#92c5de","#f4a582","#0571b0","#ca0020"))
legend("topleft", legend=c("CC", "CH", "LC", "LH"), fill = c("#92c5de","#f4a582","#0571b0","#ca0020"), bty="n")
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), pch=c(15,0,1,16), bty="n")
dev.off()

# PCo's 3 and 4
pdf(file="PCoA_ofav_host_PC4.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,3], scores[,4],col=c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))],pch=c(6,2,25,17)[as.numeric(as.factor(conditions$treat))], bg= c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))], xlab="Coordinate 3 (6.3%)", ylab="Coordinate 4 (5.3%)", main="Site")
ordiellipse(scores, conditions$site, label=F, col=c("#018571","#80cdc1","#dfc27d","#a6611a"))
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), fill = c("#018571","#80cdc1","#dfc27d","#a6611a"), bty="n")
legend("topleft", legend=c("CC", "CH", "LC", "LH"), pch=c(6,2,25,17), bty="n", pt.bg = c(NA,NA,"black",NA))
plot(scores[,3], scores[,4],col=c("#92c5de","#f4a582","#0571b0","#ca0020")[as.numeric(as.factor(conditions$treat))],pch=c(15,0,1,16)[as.numeric((as.factor(conditions$site)))], xlab="Coordinate 3 (6.3%)", ylab="Coordinate 4 (5.3%)", main="Treatment")
ordiellipse(scores, conditions$treat, label=F, col=c("#92c5de","#f4a582","#0571b0","#ca0020"))
legend("topleft", legend=c("CC", "CH", "LC", "LH"), fill = c("#92c5de","#f4a582","#0571b0","#ca0020"), bty="n")
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), pch=c(15,0,1,16), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
set.seed(12935)
ad=adonis2(t(vsd)~site*treat, data=conditions, method="manhattan", by = "terms", parallel = getOption("mc.cores"), permutations = 1e6)
ad
write.csv(ad, file = "PERMANOVA_output.csv")

# creating pie chart to represent ANOVA results
cols=c("blue","orange","grey80")
pdf(file="PERMANOVA_pie.pdf", width=6, height=6)
pie(ad$R2[1:3],labels=row.names(ad)[1:4],col=cols,main="site vs treatment")
dev.off()


#### DESEQ ####

# with multi-factor, multi-level design - using LRT
load("initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# model for the effect of treatment: (>2 factor levels => LRT)
dds$treat <- factor(dds$treat, levels = c("CC", "LC", "CH", "LH"))
dds_treat=DESeq(dds,test="LRT",reduced=~site, parallel=TRUE)

# saving all models
save(dds,dds_treat,file="realModels.RData")


#### DEGs and CONTRASTS ####

load("realModels.RData")
library(DESeq2)

# site contrasts
Rainbow_Emerald=results(dds,contrast=c("site","Rainbow","Emerald"))
summary(Rainbow_Emerald)
degs_Rainbow_Emerald=row.names(Rainbow_Emerald)[Rainbow_Emerald$padj<0.1 & !(is.na(Rainbow_Emerald$padj))]

Star_Emerald=results(dds,contrast=c("site","Star","Emerald"))
summary(Star_Emerald)
degs_Star_Emerald=row.names(Star_Emerald)[Star_Emerald$padj<0.1 & !(is.na(Star_Emerald$padj))]

MacN_Emerald=results(dds,contrast=c("site","MacN","Emerald"))
summary(MacN_Emerald)
degs_MacN_Emerald=row.names(MacN_Emerald)[MacN_Emerald$padj<0.1 & !(is.na(MacN_Emerald$padj))]

Star_Rainbow=results(dds,contrast=c("site","Star","Rainbow"))
summary(Star_Rainbow)
degs_Star_Rainbow=row.names(Star_Rainbow)[Star_Rainbow$padj<0.1 & !(is.na(Star_Rainbow$padj))]

MacN_Rainbow=results(dds,contrast=c("site","MacN","Rainbow"))
summary(MacN_Rainbow)
degs_MacN_Rainbow=row.names(MacN_Rainbow)[MacN_Rainbow$padj<0.1 & !(is.na(MacN_Rainbow$padj))]

MacN_Star=results(dds,contrast=c("site","MacN","Star"))
summary(MacN_Star)
degs_MacN_Star=row.names(MacN_Star)[MacN_Star$padj<0.1 & !(is.na(MacN_Star$padj))]

# treatment contrasts
LC_CC=results(dds,contrast=c("treat","LC","CC"))
summary(LC_CC)
degs_LC_CC=row.names(LC_CC)[LC_CC$padj<0.1 & !(is.na(LC_CC$padj))]

CH_CC=results(dds,contrast=c("treat","CH","CC"))
summary(CH_CC)
degs_CH_CC=row.names(CH_CC)[CH_CC$padj<0.1 & !(is.na(CH_CC$padj))]

LH_CC=results(dds,contrast=c("treat","LH","CC"))
summary(LH_CC)
degs_LH_CC=row.names(LH_CC)[LH_CC$padj<0.1 & !(is.na(LH_CC$padj))]

CH_LC=results(dds,contrast=c("treat","CH","LC"))
summary(CH_LC)
degs_CH_LC=row.names(CH_LC)[CH_LC$padj<0.1 & !(is.na(CH_LC$padj))]

LH_CH=results(dds,contrast=c("treat","LH","CH"))
summary(LH_CH)
degs_LH_CH=row.names(LH_CH)[LH_CH$padj<0.1 & !(is.na(LH_CH$padj))]

LH_LC=results(dds,contrast=c("treat","LH","LC"))
summary(LH_LC)
degs_LH_LC=row.names(LH_LC)[LH_LC$padj<0.1 & !(is.na(LH_LC$padj))]

save(Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star, 
     LC_CC, CH_CC, LH_CC, CH_LC, LH_CH, LH_LC, 
     degs_Rainbow_Emerald, degs_Star_Emerald, degs_MacN_Emerald, degs_Star_Rainbow, degs_MacN_Rainbow, degs_MacN_Star, 
     degs_LC_CC, degs_CH_CC, degs_LH_CC, degs_CH_LC, degs_LH_CH, degs_LH_LC, file="pvals.RData")


#### VENN DIAGRAMS ####

load("pvals.RData")
library(DESeq2)

# install.packages("VennDiagram")
library(VennDiagram)

pairwise_site=list("Star_Emerald"=degs_Star_Emerald, "MacN_Emerald"=degs_MacN_Emerald,"Star_Rainbow"=degs_Star_Rainbow, "MacN_Rainbow"=degs_MacN_Rainbow)
# pairwise comparisons among treatments
venn_site=venn.diagram(
  x = pairwise_site,
  filename=NULL,
  col = "transparent",
  fill = c("#762a83","#1b7837","#c2a5cf","#a6dba0"),
  alpha = 0.5,
  label.col = c("#9970ab","white","#5aae61","grey25","white","black","white","grey25","#762a83","white","white","white","white","#1b7837","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("#762a83","#1b7837","#9970ab","#5aae61"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_site.pdf", height=10, width=12)
grid.draw(venn_site)
dev.off()

pairwise_treat=list("LC_CC"=degs_LC_CC, "CH_CC"=degs_CH_CC, "LH_CC"=degs_LH_CC)
# pairwise comparisons among treatments
venn_treat=venn.diagram(
  x = pairwise_treat,
  filename=NULL,
  col = "transparent",
  fill = c("#0571b0","#f4a582","#ca0020"),
  alpha = 0.5,
  label.col = c("#0571b0","white","#f4a582","white","black","white","#ca0020"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("#0571b0","#f4a582","#ca0020"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0.25,0),c(0.75,0),c(0.5,1))
)
pdf(file="Venn_treatment.pdf", height=12, width=12)
grid.draw(venn_treat)
dev.off()


#### GO/KOG EXPORT ####

load("realModels.RData")
load("pvals.RData")

# site contrasts

# signed log p-values: -log(padj)* direction:
source=Rainbow_Emerald[!is.na(Rainbow_Emerald$pvalue),]
Rainbow_Emerald.p=data.frame("gene"=row.names(source))
Rainbow_Emerald.p$lpv=-log(source[,"padj"],10)
Rainbow_Emerald.p$lpv[source$stat<0]=Rainbow_Emerald.p$lpv[source$stat<0]*-1
head(Rainbow_Emerald.p)
write.csv(Rainbow_Emerald.p,file="Rainbow_Emerald_lpv.csv",row.names=F,quote=F)
save(Rainbow_Emerald.p,file="Rainbow_Emerald_lpv.RData")

source=Star_Emerald[!is.na(Star_Emerald$pvalue),]
Star_Emerald.p=data.frame("gene"=row.names(source))
Star_Emerald.p$lpv=-log(source[,"padj"],10)
Star_Emerald.p$lpv[source$stat<0]=Star_Emerald.p$lpv[source$stat<0]*-1
head(Star_Emerald.p)
write.csv(Star_Emerald.p,file="Star_Emerald_lpv.csv",row.names=F,quote=F)
save(Star_Emerald.p,file="Star_Emerald_lpv.RData")

source=MacN_Emerald[!is.na(MacN_Emerald$pvalue),]
MacN_Emerald.p=data.frame("gene"=row.names(source))
MacN_Emerald.p$lpv=-log(source[,"padj"],10)
MacN_Emerald.p$lpv[source$stat<0]=MacN_Emerald.p$lpv[source$stat<0]*-1
head(MacN_Emerald.p)
write.csv(MacN_Emerald.p,file="MacN_Emerald_lpv.csv",row.names=F,quote=F)
save(MacN_Emerald.p,file="MacN_Emerald_lpv.RData")

source=Star_Rainbow[!is.na(Star_Rainbow$pvalue),]
Star_Rainbow.p=data.frame("gene"=row.names(source))
Star_Rainbow.p$lpv=-log(source[,"padj"],10)
Star_Rainbow.p$lpv[source$stat<0]=Star_Rainbow.p$lpv[source$stat<0]*-1
head(Star_Rainbow.p)
write.csv(Star_Rainbow.p,file="Star_Rainbow_lpv.csv",row.names=F,quote=F)
save(Star_Rainbow.p,file="Star_Rainbow_lpv.RData")

source=MacN_Rainbow[!is.na(MacN_Rainbow$pvalue),]
MacN_Rainbow.p=data.frame("gene"=row.names(source))
MacN_Rainbow.p$lpv=-log(source[,"padj"],10)
MacN_Rainbow.p$lpv[source$stat<0]=MacN_Rainbow.p$lpv[source$stat<0]*-1
head(MacN_Rainbow.p)
write.csv(MacN_Rainbow.p,file="MacN_Rainbow_lpv.csv",row.names=F,quote=F)
save(MacN_Rainbow.p,file="MacN_Rainbow_lpv.RData")

source=MacN_Star[!is.na(MacN_Star$pvalue),]
MacN_Star.p=data.frame("gene"=row.names(source))
MacN_Star.p$lpv=-log(source[,"padj"],10)
MacN_Star.p$lpv[source$stat<0]=MacN_Star.p$lpv[source$stat<0]*-1
head(MacN_Star.p)
write.csv(MacN_Star.p,file="MacN_Star_lpv.csv",row.names=F,quote=F)
save(MacN_Star.p,file="MacN_Star_lpv.RData")

# treatment contrasts

# signed log p-values: -log(padj)* direction:
source=LC_CC[!is.na(LC_CC$pvalue),]
LC_CC.p=data.frame("gene"=row.names(source))
LC_CC.p$lpv=-log(source[,"padj"],10)
LC_CC.p$lpv[source$stat<0]=LC_CC.p$lpv[source$stat<0]*-1
head(LC_CC.p)
write.csv(LC_CC.p,file="LC_CC_lpv.csv",row.names=F,quote=F)
save(LC_CC.p,file="LC_CC_lpv.RData")

source=CH_CC[!is.na(CH_CC$pvalue),]
CH_CC.p=data.frame("gene"=row.names(source))
CH_CC.p$lpv=-log(source[,"padj"],10)
CH_CC.p$lpv[source$stat<0]=CH_CC.p$lpv[source$stat<0]*-1
head(CH_CC.p)
write.csv(CH_CC.p,file="CH_CC_lpv.csv",row.names=F,quote=F)
save(CH_CC.p,file="CH_CC_lpv.RData")

source=LH_CC[!is.na(LH_CC$pvalue),]
LH_CC.p=data.frame("gene"=row.names(source))
LH_CC.p$lpv=-log(source[,"padj"],10)
LH_CC.p$lpv[source$stat<0]=LH_CC.p$lpv[source$stat<0]*-1
head(LH_CC.p)
write.csv(LH_CC.p,file="LH_CC_lpv.csv",row.names=F,quote=F)
save(LH_CC.p,file="LH_CC_lpv.RData")

source=CH_LC[!is.na(CH_LC$pvalue),]
CH_LC.p=data.frame("gene"=row.names(source))
CH_LC.p$lpv=-log(source[,"padj"],10)
CH_LC.p$lpv[source$stat<0]=CH_LC.p$lpv[source$stat<0]*-1
head(CH_LC.p)
write.csv(CH_LC.p,file="CH_LC_lpv.csv",row.names=F,quote=F)
save(CH_LC.p,file="CH_LC_lpv.RData")

source=CH_LC[!is.na(CH_LC$pvalue),]
LH_LC.p=data.frame("gene"=row.names(source))
LH_LC.p$lpv=-log(source[,"padj"],10)
LH_LC.p$lpv[source$stat<0]=LH_LC.p$lpv[source$stat<0]*-1
head(LH_LC.p)
write.csv(LH_LC.p,file="LH_LC_lpv.csv",row.names=F,quote=F)
save(LH_LC.p,file="LH_LC_lpv.RData")

source=LH_CH[!is.na(LH_CH$pvalue),]
LH_CH.p=data.frame("gene"=row.names(source))
LH_CH.p$lpv=-log(source[,"padj"],10)
LH_CH.p$lpv[source$stat<0]=LH_CH.p$lpv[source$stat<0]*-1
head(LH_CH.p)
write.csv(LH_CH.p,file="LH_CH_lpv.csv",row.names=F,quote=F)
save(LH_CH.p,file="LH_CH_lpv.RData")

save(Rainbow_Emerald.p, Star_Emerald.p, MacN_Emerald.p, Star_Rainbow.p, MacN_Rainbow.p, MacN_Star.p, 
     LC_CC.p, CH_CC.p, LH_CC.p, CH_LC.p, LH_CH.p, LH_LC.p, file="exports.RData")


#### DEG MATCHING ####

library(DESeq2)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
load("exports.RData")

# This section of code does several things: 1) join -log10(pval) across treatment comparisons, 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3) adds gene annotations, and 4) then pulls on corresponding KOG classes

# Treatment
# stress treatments (LC, CH, and LH) versus control treatment (CC)
LC_CC.p %>%
  inner_join(CH_CC.p, by = c("gene" = "gene")) %>%
  inner_join(LH_CC.p, by = c("gene" = "gene")) %>%
  rename("lpv.LC_CC" = 
           lpv.x, "lpv.CH_CC" = 	
           lpv.y, "lpv.LH_CC" = 	
           lpv) %>%
  filter(abs(lpv.LC_CC) >= 1 & abs(lpv.CH_CC) >= 1 & abs(lpv.LH_CC) >= 1) %>%
  left_join(read.table(file = "../../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="stress_control", .before="gene") -> commongenes_treatment

# exporting all DEGs matching across stress vs control treatments
write.csv(commongenes_treatment, file="commongenes_treatment.csv")

# Site
# This section of code does several things: 1) join -log10(pval) across site comparisons, 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3) adds gene annotations, and 4) then pulls on corresponding KOG classes

# urban sites (Star and MacN) versus reef sites (Emerald and Rainbow)
Star_Emerald.p %>%
  inner_join(MacN_Emerald.p, by = c("gene" = "gene")) %>%
  inner_join(Star_Rainbow.p, by = c("gene" = "gene")) %>%
  inner_join(MacN_Rainbow.p, by = c("gene" = "gene")) %>%
  rename("lpv.Star_Emerald" = 
           lpv.x, "lpv.MacN_Emerald" = 	
           lpv.y, "lpv.Star_Rainbow" = 	
           lpv.x.x, "lpv.MacN_Rainbow" = 	
           lpv.y.y) %>%
  filter(abs(lpv.Star_Emerald) >= 1 & abs(lpv.MacN_Emerald) >= 1 & abs(lpv.Star_Rainbow) >= 1 & abs(lpv.MacN_Rainbow) >= 1) %>%
  left_join(read.table(file = "../../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="urban_reef", .before="gene") -> commongenes_site

# exporting all DEGs matching across stress vs control treatments
write.csv(commongenes_site, file="commongenes_site.csv")


#### HEATMAPS ####

# first removing unannotated genes
commongenes_treatment %>%
  filter(!is.na(annot)) ->  commongenes_treatment_heatmap

commongenes_site %>%
  filter(!is.na(annot)) ->  commongenes_site_heatmap

load("vsd.RData")
vsd_treatment <- subset(vsd, rownames(vsd) %in% commongenes_treatment_heatmap$gene)
vsd_site <- subset(vsd, rownames(vsd) %in% commongenes_site_heatmap$gene)

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# Treatment
# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(commongenes_treatment_heatmap$gene, commongenes_treatment_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1
pdf(file="heatmap_treatment_p0.1.pdf", height=32, width=32)
uniHeatmap(vsd=vsd_treatment,gene.names=gene_names,
           metric=-(abs(commongenes_treatment_heatmap$lpv.LH_CC)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_treatment)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 0.001
pdf(file="heatmap_treatment_p0.001.pdf", height=16, width=32)
uniHeatmap(vsd=vsd_treatment,gene.names=gene_names,
           metric=-(abs(commongenes_treatment_heatmap$lpv.LH_CC)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-3, 
           sort=c(1:ncol(vsd_treatment)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 1e-6
pdf(file="heatmap_treatment_p1e-6.pdf", height=6, width=20)
uniHeatmap(vsd=vsd_treatment,gene.names=gene_names,
           metric=-(abs(commongenes_treatment_heatmap$lpv.LH_CC)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_treatment)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# Site
# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(commongenes_site_heatmap$gene, commongenes_site_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1
pdf(file="heatmap_site_p0.1.pdf", height=20, width=24)
uniHeatmap(vsd=vsd_site,gene.names=gene_names,
           metric=-(abs(commongenes_site_heatmap$lpv.MacN_Emerald)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_site)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 0.001
pdf(file="heatmap_site_p0.001.pdf", height=12, width=24)
uniHeatmap(vsd=vsd_site,gene.names=gene_names,
           metric=-(abs(commongenes_site_heatmap$lpv.MacN_Emerald)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-3, 
           sort=c(1:ncol(vsd_site)), # overrides sorting of columns according to hierarchical clustering
            # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 1e-6
pdf(file="heatmap_site_p1e-6.pdf", height=4, width=15)
uniHeatmap(vsd=vsd_site,gene.names=gene_names,
           metric=-(abs(commongenes_site_heatmap$lpv.MacN_Emerald)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_site)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()


#### KOG MATCHING ####

# Treatment
# LC vs CC 
commongenes_treatment %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.LC_CC >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "LC_CC_up" = n) -> KOG_LC_CC_up

commongenes_treatment %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.LC_CC <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "LC_CC_down" = n) -> KOG_LC_CC_down

# CH vs CC 
commongenes_treatment %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.CH_CC >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "CH_CC_up" = n) -> KOG_CH_CC_up

commongenes_treatment %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.CH_CC <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "CH_CC_down" = n) -> KOG_CH_CC_down

# LH vs CC 
commongenes_treatment %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.LH_CC >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "LH_CC_up" = n) -> KOG_LH_CC_up

commongenes_treatment %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.LH_CC <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "LH_CC_down" = n) -> KOG_LH_CC_down

# joining all KOG class sums in a single dataframe
KOG_LC_CC_up %>%
  inner_join(KOG_CH_CC_up, by = "KOG") %>%
  inner_join(KOG_LH_CC_up, by = "KOG") %>%
  inner_join(KOG_LC_CC_down, by = "KOG") %>%
  inner_join(KOG_CH_CC_down, by = "KOG") %>%
  inner_join(KOG_LH_CC_down, by = "KOG") -> KOG_match

# melting dataframe for plotting
KOG_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_melt

# creating a custom color palette
colorCount = length(unique(KOG_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_sum <- ggplot(KOG_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Number of DEGs") +
  theme_classic()
KOG_sum
ggsave("KOG_treatment.pdf", plot= KOG_sum, width=10, height=6, units="in", dpi=300)

# Site
# Star vs Emerald 
commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.Star_Emerald >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "Star_Emerald_up" = n) -> KOG_Star_Emerald_up

commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.Star_Emerald <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "Star_Emerald_down" = n) -> KOG_Star_Emerald_down

# MacN vs Emerald 
commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.MacN_Emerald >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "MacN_Emerald_up" = n) -> KOG_MacN_Emerald_up

commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.MacN_Emerald <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "MacN_Emerald_down" = n) -> KOG_MacN_Emerald_down

# Star vs Rainbow
commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.Star_Rainbow >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "Star_Rainbow_up" = n) -> KOG_Star_Rainbow_up

commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.Star_Rainbow <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "Star_Rainbow_down" = n) -> KOG_Star_Rainbow_down

# MacN vs Rainbow
commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.MacN_Rainbow >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "MacN_Rainbow_up" = n) -> KOG_MacN_Rainbow_up

commongenes_site %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv.MacN_Rainbow <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "MacN_Rainbow_down" = n) -> KOG_MacN_Rainbow_down

# joining all KOG class sums in a single dataframe
KOG_Star_Emerald_up %>%
  inner_join(KOG_MacN_Emerald_up, by = "KOG") %>%
  inner_join(KOG_Star_Rainbow_up, by = "KOG") %>%
  inner_join(KOG_MacN_Rainbow_up, by = "KOG") %>%
  inner_join(KOG_Star_Emerald_down, by = "KOG") %>%
  inner_join(KOG_MacN_Emerald_down, by = "KOG") %>%
  inner_join(KOG_Star_Rainbow_down, by = "KOG") %>%
  inner_join(KOG_MacN_Rainbow_down, by = "KOG")-> KOG_match

# melting dataframe for plotting
KOG_match %>%
  melt(id = "KOG") %>%
  rename(comparison = variable, sum = value) -> KOG_melt

# creating a custom color palette
colorCount = length(unique(KOG_match$KOG))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))

# relative abundance plot
KOG_sum <- ggplot(KOG_melt, aes(fill = KOG, y = sum, x = comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(colorCount)) +
  labs(x = "Comparison",
       y = "Number of DEGs") +
  theme_classic()
KOG_sum
ggsave("KOG_site.pdf", plot= KOG_sum, width=14, height=6, units="in", dpi=300)


#### CHERRY PICKING ####

library(DESeq2)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
load("exports.RData")

MacN_Emerald.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'heat shock|carbonic|digest|pattern recognition|respiration|immun|NF-kappaB|TGF-beta|peroxidas|protein tyrosine kinase|WD repeat-containing protein|fibrinogen|apoptosis|stress|extracellular matrix')) -> cherrypicking_site
write.csv(cherrypicking_site, file = "cherrypicking_MacN_Emerald.csv")

LH_CC.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'heat shock|carbonic|digest|pattern recognition|respiration|immun|NF-kappaB|TGF-beta|peroxidas|protein tyrosine kinase|WD repeat-containing protein|fibrinogen|apoptosis|stress|extracellular matrix')) -> cherrypicking_treat
write.csv(cherrypicking_treat, file = "cherrypicking_LH_CC.csv")


#### CHERRY SITE EXPORTS ####

library(DESeq2)
library(ggpubr)
load("realModels.RData")

# exporting counts of specific genes from cherry picking
Ofaveolata001579 <- plotCounts(dds, gene="Ofaveolata001579", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata001579, file = "Ofaveolata001579_site.csv")

Ofaveolata003869 <- plotCounts(dds, gene="Ofaveolata003869", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata003869, file = "Ofaveolata003869_site.csv")

Ofaveolata003873 <- plotCounts(dds, gene="Ofaveolata003873", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata003873, file = "Ofaveolata003873_site.csv")

Ofaveolata004069 <- plotCounts(dds, gene="Ofaveolata004069", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata004069, file = "Ofaveolata004069_site.csv")

Ofaveolata004373 <- plotCounts(dds, gene="Ofaveolata004373", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata004373, file = "Ofaveolata004373_site.csv")

Ofaveolata006068 <- plotCounts(dds, gene="Ofaveolata006068", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata006068, file = "Ofaveolata006068_site.csv")

Ofaveolata007015 <- plotCounts(dds, gene="Ofaveolata007015", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata007015, file = "Ofaveolata007015_site.csv")

Ofaveolata007419 <- plotCounts(dds, gene="Ofaveolata007419", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata007419, file = "Ofaveolata007419_site.csv")

Ofaveolata007441 <- plotCounts(dds, gene="Ofaveolata007441", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata007441, file = "Ofaveolata007441_site.csv")

Ofaveolata007495 <- plotCounts(dds, gene="Ofaveolata007495", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata007495, file = "Ofaveolata007495_site.csv")

Ofaveolata007753 <- plotCounts(dds, gene="Ofaveolata007753", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata007753, file = "Ofaveolata007753_site.csv")

Ofaveolata007833 <- plotCounts(dds, gene="Ofaveolata007833", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata007833, file = "Ofaveolata007833_site.csv")

Ofaveolata008267 <- plotCounts(dds, gene="Ofaveolata008267", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata008267, file = "Ofaveolata008267_site.csv")

Ofaveolata008270 <- plotCounts(dds, gene="Ofaveolata008270", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata008270, file = "Ofaveolata008270_site.csv")

Ofaveolata008603 <- plotCounts(dds, gene="Ofaveolata008603", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata008603, file = "Ofaveolata008603_site.csv")

Ofaveolata009657 <- plotCounts(dds, gene="Ofaveolata009657", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata009657, file = "Ofaveolata009657_site.csv")

Ofaveolata009704 <- plotCounts(dds, gene="Ofaveolata009704", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata009704, file = "Ofaveolata009704_site.csv")

Ofaveolata015105 <- plotCounts(dds, gene="Ofaveolata015105", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata015105, file = "Ofaveolata015105_site.csv")

Ofaveolata015374 <- plotCounts(dds, gene="Ofaveolata015374", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata015374, file = "Ofaveolata015374_site.csv")

Ofaveolata016018 <- plotCounts(dds, gene="Ofaveolata016018", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata016018, file = "Ofaveolata016018_site.csv")

Ofaveolata016019 <- plotCounts(dds, gene="Ofaveolata016019", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata016019, file = "Ofaveolata016019_site.csv")

Ofaveolata018959 <- plotCounts(dds, gene="Ofaveolata018959", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata018959, file = "Ofaveolata018959_site.csv")

Ofaveolata019263 <- plotCounts(dds, gene="Ofaveolata019263", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata019263, file = "Ofaveolata019263_site.csv")

Ofaveolata019402 <- plotCounts(dds, gene="Ofaveolata019402", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata019402, file = "Ofaveolata019402_site.csv")

Ofaveolata020014 <- plotCounts(dds, gene="Ofaveolata020014", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata020014, file = "Ofaveolata020014_site.csv")

Ofaveolata021096 <- plotCounts(dds, gene="Ofaveolata021096", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata021096, file = "Ofaveolata021096_site.csv")

Ofaveolata021961 <- plotCounts(dds, gene="Ofaveolata021961", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata021961, file = "Ofaveolata021961_site.csv")

Ofaveolata022243 <- plotCounts(dds, gene="Ofaveolata022243", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata022243, file = "Ofaveolata022243_site.csv")

Ofaveolata022997 <- plotCounts(dds, gene="Ofaveolata022997", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata022997, file = "Ofaveolata022997_site.csv")

Ofaveolata024007 <- plotCounts(dds, gene="Ofaveolata024007", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata024007, file = "Ofaveolata024007_site.csv")

Ofaveolata024322 <- plotCounts(dds, gene="Ofaveolata024322", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata024322, file = "Ofaveolata024322_site.csv")

Ofaveolata024334 <- plotCounts(dds, gene="Ofaveolata024334", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata024334, file = "Ofaveolata024334_site.csv")

Ofaveolata024336 <- plotCounts(dds, gene="Ofaveolata024336", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata024336, file = "Ofaveolata024336_site.csv")

Ofaveolata024338 <- plotCounts(dds, gene="Ofaveolata024338", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata024338, file = "Ofaveolata024338_site.csv")

Ofaveolata025025 <- plotCounts(dds, gene="Ofaveolata025025", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata025025, file = "Ofaveolata025025_site.csv")

Ofaveolata025300 <- plotCounts(dds, gene="Ofaveolata025300", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata025300, file = "Ofaveolata025300_site.csv")

Ofaveolata025301 <- plotCounts(dds, gene="Ofaveolata025301", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata025301, file = "Ofaveolata025301_site.csv")

Ofaveolata026803 <- plotCounts(dds, gene="Ofaveolata026803", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata026803, file = "Ofaveolata026803_site.csv")

Ofaveolata027116 <- plotCounts(dds, gene="Ofaveolata027116", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata027116, file = "Ofaveolata027116_site.csv")

Ofaveolata027542 <- plotCounts(dds, gene="Ofaveolata027542", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata027542, file = "Ofaveolata027542_site.csv")

Ofaveolata028116 <- plotCounts(dds, gene="Ofaveolata028116", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata028116, file = "Ofaveolata028116_site.csv")

Ofaveolata030665 <- plotCounts(dds, gene="Ofaveolata030665", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata030665, file = "Ofaveolata030665_site.csv")

Ofaveolata030749 <- plotCounts(dds, gene="Ofaveolata030749", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata030749, file = "Ofaveolata030749_site.csv")

Ofaveolata030946 <- plotCounts(dds, gene="Ofaveolata030946", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata030946, file = "Ofaveolata030946_site.csv")

Ofaveolata031073 <- plotCounts(dds, gene="Ofaveolata031073", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata031073, file = "Ofaveolata031073_site.csv")

Ofaveolata032411 <- plotCounts(dds, gene="Ofaveolata032411", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata032411, file = "Ofaveolata032411_site.csv")

Ofaveolata033022 <- plotCounts(dds, gene="Ofaveolata033022", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata033022, file = "Ofaveolata033022_site.csv")

Ofaveolata033029 <- plotCounts(dds, gene="Ofaveolata033029", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata033029, file = "Ofaveolata033029_site.csv")

Ofaveolata033030 <- plotCounts(dds, gene="Ofaveolata033030", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata033030, file = "Ofaveolata033030_site.csv")

Ofaveolata033153 <- plotCounts(dds, gene="Ofaveolata033153", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata033153, file = "Ofaveolata033153_site.csv")

Ofaveolata035052 <- plotCounts(dds, gene="Ofaveolata035052", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata035052, file = "Ofaveolata035052_site.csv")

Ofaveolata035203 <- plotCounts(dds, gene="Ofaveolata035203", intgroup="site", returnData=TRUE)
write.csv(Ofaveolata035203, file = "Ofaveolata035203_site.csv")


#### CHERRY TREATMENT EXPORTS ####

library(DESeq2)
library(ggpubr)
load("realModels.RData")

# exporting counts of specific genes from cherry picking
Ofaveolata001477 <- plotCounts(dds, gene="Ofaveolata001477", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata001477, file = "Ofaveolata001477_treat.csv")

Ofaveolata005172 <- plotCounts(dds, gene="Ofaveolata005172", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata005172, file = "Ofaveolata005172_treat.csv")

Ofaveolata007015 <- plotCounts(dds, gene="Ofaveolata007015", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata007015, file = "Ofaveolata007015_treat.csv")

Ofaveolata007495 <- plotCounts(dds, gene="Ofaveolata007495", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata007495, file = "Ofaveolata007495_treat.csv")

Ofaveolata007753 <- plotCounts(dds, gene="Ofaveolata007753", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata007753, file = "Ofaveolata007753_treat.csv")

Ofaveolata007848 <- plotCounts(dds, gene="Ofaveolata007848", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata007848, file = "Ofaveolata007848_treat.csv")

Ofaveolata008267 <- plotCounts(dds, gene="Ofaveolata008267", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata008267, file = "Ofaveolata008267_treat.csv")

Ofaveolata008270 <- plotCounts(dds, gene="Ofaveolata008270", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata008270, file = "Ofaveolata008270_treat.csv")

Ofaveolata010796 <- plotCounts(dds, gene="Ofaveolata010796", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata010796, file = "Ofaveolata010796_treat.csv")

Ofaveolata011375 <- plotCounts(dds, gene="Ofaveolata011375", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata011375, file = "Ofaveolata011375_treat.csv")

Ofaveolata013074 <- plotCounts(dds, gene="Ofaveolata013074", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata013074, file = "Ofaveolata013074_treat.csv")

Ofaveolata016018 <- plotCounts(dds, gene="Ofaveolata016018", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata016018, file = "Ofaveolata016018_treat.csv")

Ofaveolata016019 <- plotCounts(dds, gene="Ofaveolata016019", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata016019, file = "Ofaveolata016019_treat.csv")

Ofaveolata016020 <- plotCounts(dds, gene="Ofaveolata016020", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata016020, file = "Ofaveolata016020_treat.csv")

Ofaveolata016021 <- plotCounts(dds, gene="Ofaveolata016021", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata016021, file = "Ofaveolata016021_treat.csv")

Ofaveolata018959 <- plotCounts(dds, gene="Ofaveolata018959", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata018959, file = "Ofaveolata018959_treat.csv")

Ofaveolata019944 <- plotCounts(dds, gene="Ofaveolata019944", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata019944, file = "Ofaveolata019944_treat.csv")

Ofaveolata021361 <- plotCounts(dds, gene="Ofaveolata021361", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata021361, file = "Ofaveolata021361_treat.csv")

Ofaveolata022997 <- plotCounts(dds, gene="Ofaveolata022997", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata022997, file = "Ofaveolata022997_treat.csv")

Ofaveolata024321 <- plotCounts(dds, gene="Ofaveolata024321", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024321, file = "Ofaveolata024321_treat.csv")

Ofaveolata024323 <- plotCounts(dds, gene="Ofaveolata024323", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024323, file = "Ofaveolata024323_treat.csv")

Ofaveolata024326 <- plotCounts(dds, gene="Ofaveolata024326", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024326, file = "Ofaveolata024326_treat.csv")

Ofaveolata024334 <- plotCounts(dds, gene="Ofaveolata024334", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024334, file = "Ofaveolata024334_treat.csv")

Ofaveolata024336 <- plotCounts(dds, gene="Ofaveolata024336", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024336, file = "Ofaveolata024336_treat.csv")

Ofaveolata024338 <- plotCounts(dds, gene="Ofaveolata024338", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024338, file = "Ofaveolata024338_treat.csv")

Ofaveolata024342 <- plotCounts(dds, gene="Ofaveolata024342", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata024342, file = "Ofaveolata024342_treat.csv")

Ofaveolata025300 <- plotCounts(dds, gene="Ofaveolata025300", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata025300, file = "Ofaveolata025300_treat.csv")

Ofaveolata026158 <- plotCounts(dds, gene="Ofaveolata026158", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata026158, file = "Ofaveolata026158_treat.csv")

Ofaveolata027268 <- plotCounts(dds, gene="Ofaveolata027268", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata027268, file = "Ofaveolata027268_treat.csv")

Ofaveolata027534 <- plotCounts(dds, gene="Ofaveolata027534", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata027534, file = "Ofaveolata027534_treat.csv")

Ofaveolata027542 <- plotCounts(dds, gene="Ofaveolata027542", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata027542, file = "Ofaveolata027542_treat.csv")

Ofaveolata027545 <- plotCounts(dds, gene="Ofaveolata027545", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata027545, file = "Ofaveolata027545_treat.csv")

Ofaveolata028699 <- plotCounts(dds, gene="Ofaveolata028699", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata028699, file = "Ofaveolata028699_treat.csv")

Ofaveolata032821 <- plotCounts(dds, gene="Ofaveolata032821", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata032821, file = "Ofaveolata032821_treat.csv")

Ofaveolata035792 <- plotCounts(dds, gene="Ofaveolata035792", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata035792, file = "Ofaveolata035792_treat.csv")

Ofaveolata035793 <- plotCounts(dds, gene="Ofaveolata035793", intgroup="treat", returnData=TRUE)
write.csv(Ofaveolata035793, file = "Ofaveolata035793_treat.csv")
