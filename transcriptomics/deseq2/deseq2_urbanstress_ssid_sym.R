#### PACKAGES ####

# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
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
library(dplyr)

#read in counts
counts = read.table("../../../../raw/ssid/allcounts_sym.txt")

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
design = read.csv("../../../../raw/design_ssid.csv", head=TRUE)
design
design$site <- as.factor(design$site)
design$site <- as.factor(design$site)
design$replicate <- as.factor(design$replicate)
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
outs=c(7,8,26,27,41,46,47,66,67)
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
pdf(file="heatmap_ssid_sym.pdf", width=15, height=15)
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
pdf(file="PCoA_ssid_sym.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))],pch=c(6,2,25,17)[as.numeric(as.factor(conditions$treat))], bg= c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1 (37.0%)", ylab="Coordinate 2 (4.9%)", main="Site")
ordiellipse(scores, conditions$site, label=F, col=c("#018571","#80cdc1","#dfc27d","#a6611a"))
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), fill = c("#018571","#80cdc1","#dfc27d","#a6611a"), bty="n")
legend("topleft", legend=c("CC", "CH", "LC", "LH"), pch=c(6,2,25,17), bty="n", pt.bg = c(NA,NA,"black",NA))
plot(scores[,1], scores[,2],col=c("#92c5de","#f4a582","#0571b0","#ca0020")[as.numeric(as.factor(conditions$treat))],pch=c(15,0,1,16)[as.numeric((as.factor(conditions$site)))], xlab="Coordinate 1 (37.0%)", ylab="Coordinate 2 (4.9%)", main="Treatment")
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
set.seed(53284)
ad=adonis2(t(vsd)~site*treat, data=conditions, method="manhattan", permutations=1e6)
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
Rainbow_Emerald.p=data.frame("gene"=row.names(source))
Rainbow_Emerald.p$lpv=-log(source[,"padj"],10)
Rainbow_Emerald.p$lpv[source$stat<0]=Rainbow_Emerald.p$lpv[source$stat<0]*-1
head(Rainbow_Emerald.p)
write.csv(Rainbow_Emerald.p,file="Rainbow_Emerald_lpv.csv",row.names=F,quote=F)
save(Rainbow_Emerald.p,file="Rainbow_Emerald_lpv.RData")

Star_Emerald.p=data.frame("gene"=row.names(source))
Star_Emerald.p$lpv=-log(source[,"padj"],10)
Star_Emerald.p$lpv[source$stat<0]=Star_Emerald.p$lpv[source$stat<0]*-1
head(Star_Emerald.p)
write.csv(Star_Emerald.p,file="Star_Emerald_lpv.csv",row.names=F,quote=F)
save(Star_Emerald.p,file="Star_Emerald_lpv.RData")

MacN_Emerald.p=data.frame("gene"=row.names(source))
MacN_Emerald.p$lpv=-log(source[,"padj"],10)
MacN_Emerald.p$lpv[source$stat<0]=MacN_Emerald.p$lpv[source$stat<0]*-1
head(MacN_Emerald.p)
write.csv(MacN_Emerald.p,file="MacN_Emerald_lpv.csv",row.names=F,quote=F)
save(MacN_Emerald.p,file="MacN_Emerald_lpv.RData")

Star_Rainbow.p=data.frame("gene"=row.names(source))
Star_Rainbow.p$lpv=-log(source[,"padj"],10)
Star_Rainbow.p$lpv[source$stat<0]=Star_Rainbow.p$lpv[source$stat<0]*-1
head(Star_Rainbow.p)
write.csv(Star_Rainbow.p,file="Star_Rainbow_lpv.csv",row.names=F,quote=F)
save(Star_Rainbow.p,file="Star_Rainbow_lpv.RData")

MacN_Rainbow.p=data.frame("gene"=row.names(source))
MacN_Rainbow.p$lpv=-log(source[,"padj"],10)
MacN_Rainbow.p$lpv[source$stat<0]=MacN_Rainbow.p$lpv[source$stat<0]*-1
head(MacN_Rainbow.p)
write.csv(MacN_Rainbow.p,file="MacN_Rainbow_lpv.csv",row.names=F,quote=F)
save(MacN_Rainbow.p,file="MacN_Rainbow_lpv.RData")

MacN_Star.p=data.frame("gene"=row.names(source))
MacN_Star.p$lpv=-log(source[,"padj"],10)
MacN_Star.p$lpv[source$stat<0]=MacN_Star.p$lpv[source$stat<0]*-1
head(MacN_Star.p)
write.csv(MacN_Star.p,file="MacN_Star_lpv.csv",row.names=F,quote=F)
save(MacN_Star.p,file="MacN_Star_lpv.RData")

# treatment contrasts

# signed log p-values: -log(padj)* direction:
LC_CC.p=data.frame("gene"=row.names(source))
LC_CC.p$lpv=-log(source[,"padj"],10)
LC_CC.p$lpv[source$stat<0]=LC_CC.p$lpv[source$stat<0]*-1
head(LC_CC.p)
write.csv(LC_CC.p,file="LC_CC_lpv.csv",row.names=F,quote=F)
save(LC_CC.p,file="LC_CC_lpv.RData")

CH_CC.p=data.frame("gene"=row.names(source))
CH_CC.p$lpv=-log(source[,"padj"],10)
CH_CC.p$lpv[source$stat<0]=CH_CC.p$lpv[source$stat<0]*-1
head(CH_CC.p)
write.csv(CH_CC.p,file="CH_CC_lpv.csv",row.names=F,quote=F)
save(CH_CC.p,file="CH_CC_lpv.RData")

LH_CC.p=data.frame("gene"=row.names(source))
LH_CC.p$lpv=-log(source[,"padj"],10)
LH_CC.p$lpv[source$stat<0]=LH_CC.p$lpv[source$stat<0]*-1
head(LH_CC.p)
write.csv(LH_CC.p,file="LH_CC_lpv.csv",row.names=F,quote=F)
save(LH_CC.p,file="LH_CC_lpv.RData")

CH_LC.p=data.frame("gene"=row.names(source))
CH_LC.p$lpv=-log(source[,"padj"],10)
CH_LC.p$lpv[source$stat<0]=CH_LC.p$lpv[source$stat<0]*-1
head(CH_LC.p)
write.csv(CH_LC.p,file="CH_LC_lpv.csv",row.names=F,quote=F)
save(CH_LC.p,file="CH_LC_lpv.RData")

LH_LC.p=data.frame("gene"=row.names(source))
LH_LC.p$lpv=-log(source[,"padj"],10)
LH_LC.p$lpv[source$stat<0]=LH_LC.p$lpv[source$stat<0]*-1
head(LH_LC.p)
write.csv(LH_LC.p,file="LH_LC_lpv.csv",row.names=F,quote=F)
save(LH_LC.p,file="LH_LC_lpv.RData")

LH_CH.p=data.frame("gene"=row.names(source))
LH_CH.p$lpv=-log(source[,"padj"],10)
LH_CH.p$lpv[source$stat<0]=LH_CH.p$lpv[source$stat<0]*-1
head(LH_CH.p)
write.csv(LH_CH.p,file="LH_CH_lpv.csv",row.names=F,quote=F)
save(LH_CH.p,file="LH_CH_lpv.RData")



