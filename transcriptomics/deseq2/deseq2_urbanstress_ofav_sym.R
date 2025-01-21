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
library(dplyr)

#read in counts
counts = read.table("../../../../../raw/ofav/orthogroup/allcounts_sym.txt")

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
outs=c(22,56)
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
pdf(file="heatmap_ofav_sym.pdf", width=15, height=15)
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
pdf(file="PCoA_ofav_sym.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))],pch=c(6,2,25,17)[as.numeric(as.factor(conditions$treat))], bg= c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1 (69.3%)", ylab="Coordinate 2 (5.5%)", main="Site")
ordiellipse(scores, conditions$site, label=F, col=c("#018571","#80cdc1","#dfc27d","#a6611a"))
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), fill = c("#018571","#80cdc1","#dfc27d","#a6611a"), bty="n")
legend("topleft", legend=c("CC", "CH", "LC", "LH"), pch=c(6,2,25,17), bty="n", pt.bg = c(NA,NA,"black",NA))
plot(scores[,1], scores[,2],col=c("#92c5de","#f4a582","#0571b0","#ca0020")[as.numeric(as.factor(conditions$treat))],pch=c(15,0,1,16)[as.numeric((as.factor(conditions$site)))], xlab="Coordinate 1 (69.3%)", ylab="Coordinate 2 (5.5%)", main="Treatment")
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
set.seed(25786)
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
library(dplyr)
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
  left_join(read.table(file = "../../../../../Annotations/symD/shoguchi/pcli_host/Durusdinium_ortho2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(ortho = V1,
                     annot = V2) %>%
              distinct(ortho, .keep_all = TRUE) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "ortho")) %>%
  left_join(read.table(file = "../../../../../Annotations/symD/shoguchi/pcli_host/Durusdinium_ortho2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(ortho = V1,
                     KOG = V2) %>%
              distinct(ortho, .keep_all = TRUE) %>%
              dplyr::select(-V1, -V2),
            by = c("gene" = "ortho")) %>%
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
  left_join(read.table(file = "../../../../../Annotations/symD/shoguchi/pcli_host/Durusdinium_ortho2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(ortho = V1,
                     annot = V2) %>%
              distinct(ortho, .keep_all = TRUE) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "ortho")) %>%
  left_join(read.table(file = "../../../../../Annotations/symD/shoguchi/pcli_host/Durusdinium_ortho2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(ortho = V1,
                     KOG = V2) %>%
              distinct(ortho, .keep_all = TRUE) %>%
              dplyr::select(-V1, -V2),
              by = c("gene" = "ortho")) %>%
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
pdf(file="heatmap_treatment_p0.1.pdf", height=1.5, width=18)
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

# Site
# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(commongenes_site_heatmap$gene, commongenes_site_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1
pdf(file="heatmap_site_p0.1.pdf", height=96, width=32)
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

# p < 1e-6
pdf(file="heatmap_site_p1e-6.pdf", height=12, width=32)
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
  inner_join(KOG_LH_CC_up, by = "KOG") -> KOG_match
  # inner_join(KOG_LC_CC_down, by = "KOG") %>%
  # inner_join(KOG_CH_CC_down, by = "KOG") %>%
  # inner_join(KOG_LH_CC_down, by = "KOG") 

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
ggsave("KOG_treatment.pdf", plot= KOG_sum, width=7, height=6, units="in", dpi=300)

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
  left_join(read.table(file = "../../../../../Annotations/orthogroup/Durusdinium_ortho2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'heat shock|carbonic|digest|pattern recognition|respiration|immun|NF-kappaB|TGF-beta|peroxidas|protein tyrosine kinase|WD repeat-containing protein|fibrinogen|apoptosis|stress|extracellular matrix|photsynthe|thylakoid')) -> cherrypicking_site
write.csv(cherrypicking_site, file = "cherrypicking_MacN_Emerald.csv")

LH_CC.p %>%
  filter(abs(lpv) >= 1) %>%
  left_join(read.table(file = "../../../../../Annotations/orthogroup/Durusdinium_ortho2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  filter(str_detect(annot, 'heat shock|carbonic|digest|pattern recognition|respiration|immun|NF-kappaB|TGF-beta|peroxidas|protein tyrosine kinase|WD repeat-containing protein|fibrinogen|apoptosis|stress|extracellular matrix|photsynthe|thylakoid')) -> cherrypicking_treat
write.csv(cherrypicking_treat, file = "cherrypicking_LH_CC.csv")


#### CHERRY SITE EXPORTS ####

library(DESeq2)
library(ggpubr)
load("realModels.RData")

# exporting counts of specific genes from cherry picking
Durusdinium11872 <- plotCounts(dds, gene="Durusdinium11872", intgroup="site", returnData=TRUE)
write.csv(Durusdinium11872, file = "Durusdinium11872_site.csv")

Durusdinium11932 <- plotCounts(dds, gene="Durusdinium11932", intgroup="site", returnData=TRUE)
write.csv(Durusdinium11932, file = "Durusdinium11932_site.csv")

Durusdinium12554 <- plotCounts(dds, gene="Durusdinium12554", intgroup="site", returnData=TRUE)
write.csv(Durusdinium12554, file = "Durusdinium12554_site.csv")

Durusdinium13344 <- plotCounts(dds, gene="Durusdinium13344", intgroup="site", returnData=TRUE)
write.csv(Durusdinium13344, file = "Durusdinium13344_site.csv")

Durusdinium16441 <- plotCounts(dds, gene="Durusdinium16441", intgroup="site", returnData=TRUE)
write.csv(Durusdinium16441, file = "Durusdinium16441_site.csv")

Durusdinium17451 <- plotCounts(dds, gene="Durusdinium17451", intgroup="site", returnData=TRUE)
write.csv(Durusdinium17451, file = "Durusdinium17451_site.csv")

Durusdinium17786 <- plotCounts(dds, gene="Durusdinium17786", intgroup="site", returnData=TRUE)
write.csv(Durusdinium17786, file = "Durusdinium17786_site.csv")

Durusdinium18167 <- plotCounts(dds, gene="Durusdinium18167", intgroup="site", returnData=TRUE)
write.csv(Durusdinium18167, file = "Durusdinium18167_site.csv")

Durusdinium19864 <- plotCounts(dds, gene="Durusdinium19864", intgroup="site", returnData=TRUE)
write.csv(Durusdinium19864, file = "Durusdinium19864_site.csv")

Durusdinium19903 <- plotCounts(dds, gene="Durusdinium19903", intgroup="site", returnData=TRUE)
write.csv(Durusdinium19903, file = "Durusdinium19903_site.csv")

Durusdinium21945 <- plotCounts(dds, gene="Durusdinium21945", intgroup="site", returnData=TRUE)
write.csv(Durusdinium21945, file = "Durusdinium21945_site.csv")

Durusdinium22479 <- plotCounts(dds, gene="Durusdinium22479", intgroup="site", returnData=TRUE)
write.csv(Durusdinium22479, file = "Durusdinium22479_site.csv")

Durusdinium23177 <- plotCounts(dds, gene="Durusdinium23177", intgroup="site", returnData=TRUE)
write.csv(Durusdinium23177, file = "Durusdinium23177_site.csv")

Durusdinium23955 <- plotCounts(dds, gene="Durusdinium23955", intgroup="site", returnData=TRUE)
write.csv(Durusdinium23955, file = "Durusdinium23955_site.csv")

Durusdinium25695 <- plotCounts(dds, gene="Durusdinium25695", intgroup="site", returnData=TRUE)
write.csv(Durusdinium25695, file = "Durusdinium25695_site.csv")

Durusdinium26534 <- plotCounts(dds, gene="Durusdinium26534", intgroup="site", returnData=TRUE)
write.csv(Durusdinium26534, file = "Durusdinium26534_site.csv")

Durusdinium26834 <- plotCounts(dds, gene="Durusdinium26834", intgroup="site", returnData=TRUE)
write.csv(Durusdinium26834, file = "Durusdinium26834_site.csv")

Durusdinium27173 <- plotCounts(dds, gene="Durusdinium27173", intgroup="site", returnData=TRUE)
write.csv(Durusdinium27173, file = "Durusdinium27173_site.csv")

Durusdinium27547 <- plotCounts(dds, gene="Durusdinium27547", intgroup="site", returnData=TRUE)
write.csv(Durusdinium27547, file = "Durusdinium27547_site.csv")

Durusdinium28251 <- plotCounts(dds, gene="Durusdinium28251", intgroup="site", returnData=TRUE)
write.csv(Durusdinium28251, file = "Durusdinium28251_site.csv")

Durusdinium32490 <- plotCounts(dds, gene="Durusdinium32490", intgroup="site", returnData=TRUE)
write.csv(Durusdinium32490, file = "Durusdinium32490_site.csv")

Durusdinium32556 <- plotCounts(dds, gene="Durusdinium32556", intgroup="site", returnData=TRUE)
write.csv(Durusdinium32556, file = "Durusdinium32556_site.csv")

Durusdinium33723 <- plotCounts(dds, gene="Durusdinium33723", intgroup="site", returnData=TRUE)
write.csv(Durusdinium33723, file = "Durusdinium33723_site.csv")

Durusdinium34070 <- plotCounts(dds, gene="Durusdinium34070", intgroup="site", returnData=TRUE)
write.csv(Durusdinium34070, file = "Durusdinium34070_site.csv")

Durusdinium34946 <- plotCounts(dds, gene="Durusdinium34946", intgroup="site", returnData=TRUE)
write.csv(Durusdinium34946, file = "Durusdinium34946_site.csv")

Durusdinium35558 <- plotCounts(dds, gene="Durusdinium35558", intgroup="site", returnData=TRUE)
write.csv(Durusdinium35558, file = "Durusdinium35558_site.csv")

Durusdinium36943 <- plotCounts(dds, gene="Durusdinium36943", intgroup="site", returnData=TRUE)
write.csv(Durusdinium36943, file = "Durusdinium36943_site.csv")

Durusdinium37005 <- plotCounts(dds, gene="Durusdinium37005", intgroup="site", returnData=TRUE)
write.csv(Durusdinium37005, file = "Durusdinium37005_site.csv")

Durusdinium37662 <- plotCounts(dds, gene="Durusdinium37662", intgroup="site", returnData=TRUE)
write.csv(Durusdinium37662, file = "Durusdinium37662_site.csv")

Durusdinium38393 <- plotCounts(dds, gene="Durusdinium38393", intgroup="site", returnData=TRUE)
write.csv(Durusdinium38393, file = "Durusdinium38393_site.csv")

Durusdinium38691 <- plotCounts(dds, gene="Durusdinium38691", intgroup="site", returnData=TRUE)
write.csv(Durusdinium38691, file = "Durusdinium38691_site.csv")

Durusdinium40014 <- plotCounts(dds, gene="Durusdinium40014", intgroup="site", returnData=TRUE)
write.csv(Durusdinium40014, file = "Durusdinium40014_site.csv")

Durusdinium40879 <- plotCounts(dds, gene="Durusdinium40879", intgroup="site", returnData=TRUE)
write.csv(Durusdinium40879, file = "Durusdinium40879_site.csv")

Durusdinium41126 <- plotCounts(dds, gene="Durusdinium41126", intgroup="site", returnData=TRUE)
write.csv(Durusdinium41126, file = "Durusdinium41126_site.csv")

Durusdinium42291 <- plotCounts(dds, gene="Durusdinium42291", intgroup="site", returnData=TRUE)
write.csv(Durusdinium42291, file = "Durusdinium42291_site.csv")

Durusdinium43213 <- plotCounts(dds, gene="Durusdinium43213", intgroup="site", returnData=TRUE)
write.csv(Durusdinium43213, file = "Durusdinium43213_site.csv")

Durusdinium44555 <- plotCounts(dds, gene="Durusdinium44555", intgroup="site", returnData=TRUE)
write.csv(Durusdinium44555, file = "Durusdinium44555_site.csv")

Durusdinium45079 <- plotCounts(dds, gene="Durusdinium45079", intgroup="site", returnData=TRUE)
write.csv(Durusdinium45079, file = "Durusdinium45079_site.csv")

Durusdinium4559 <- plotCounts(dds, gene="Durusdinium4559", intgroup="site", returnData=TRUE)
write.csv(Durusdinium4559, file = "Durusdinium4559_site.csv")

Durusdinium46326 <- plotCounts(dds, gene="Durusdinium46326", intgroup="site", returnData=TRUE)
write.csv(Durusdinium46326, file = "Durusdinium46326_site.csv")

Durusdinium5039 <- plotCounts(dds, gene="Durusdinium5039", intgroup="site", returnData=TRUE)
write.csv(Durusdinium5039, file = "Durusdinium5039_site.csv")

Durusdinium50438 <- plotCounts(dds, gene="Durusdinium50438", intgroup="site", returnData=TRUE)
write.csv(Durusdinium50438, file = "Durusdinium50438_site.csv")

Durusdinium52737 <- plotCounts(dds, gene="Durusdinium52737", intgroup="site", returnData=TRUE)
write.csv(Durusdinium52737, file = "Durusdinium52737_site.csv")

Durusdinium54311 <- plotCounts(dds, gene="Durusdinium54311", intgroup="site", returnData=TRUE)
write.csv(Durusdinium54311, file = "Durusdinium54311_site.csv")

Durusdinium7837 <- plotCounts(dds, gene="Durusdinium7837", intgroup="site", returnData=TRUE)
write.csv(Durusdinium7837, file = "Durusdinium7837_site.csv")


#### CHERRY TREATMENT EXPORTS ####

library(DESeq2)
library(ggpubr)
load("realModels.RData")

# exporting counts of specific genes from cherry picking
Durusdinium16264 <- plotCounts(dds, gene="Durusdinium16264", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium16264, file = "Durusdinium16264_treat.csv")

Durusdinium17786 <- plotCounts(dds, gene="Durusdinium17786", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium17786, file = "Durusdinium17786_treat.csv")

Durusdinium18167 <- plotCounts(dds, gene="Durusdinium18167", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium18167, file = "Durusdinium18167_treat.csv")

Durusdinium19864 <- plotCounts(dds, gene="Durusdinium19864", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium19864, file = "Durusdinium19864_treat.csv")

Durusdinium28251 <- plotCounts(dds, gene="Durusdinium28251", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium28251, file = "Durusdinium28251_treat.csv")

Durusdinium32490 <- plotCounts(dds, gene="Durusdinium32490", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium32490, file = "Durusdinium32490_treat.csv")

Durusdinium32556 <- plotCounts(dds, gene="Durusdinium32556", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium32556, file = "Durusdinium32556_treat.csv")

Durusdinium33723 <- plotCounts(dds, gene="Durusdinium33723", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium33723, file = "Durusdinium33723_treat.csv")

Durusdinium37662 <- plotCounts(dds, gene="Durusdinium37662", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium37662, file = "Durusdinium37662_treat.csv")

Durusdinium40014 <- plotCounts(dds, gene="Durusdinium40014", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium40014, file = "Durusdinium40014_treat.csv")

Durusdinium45079 <- plotCounts(dds, gene="Durusdinium45079", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium45079, file = "Durusdinium45079_treat.csv")

Durusdinium498 <- plotCounts(dds, gene="Durusdinium498", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium498, file = "Durusdinium498_treat.csv")

Durusdinium532 <- plotCounts(dds, gene="Durusdinium532", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium532, file = "Durusdinium532_treat.csv")

Durusdinium54280 <- plotCounts(dds, gene="Durusdinium54280", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium54280, file = "Durusdinium54280_treat.csv")

Durusdinium54311 <- plotCounts(dds, gene="Durusdinium54311", intgroup="treat", returnData=TRUE)
write.csv(Durusdinium54311, file = "Durusdinium54311_treat.csv")
