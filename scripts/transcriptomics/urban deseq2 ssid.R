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
counts = read.table("../../data/transcriptomics/ssid/allcounts_host.txt")

# how many genes we have total?
nrow(counts) 
ncol(counts)

# how does the data look? 
head(counts)

keep <- rowSums(counts) >= 10
countData <- counts[keep,]
nrow(countData)
ncol(countData)
write.csv(countData, file="../../outputs/transcriptomics/ssid/deseq2/countData.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = counts[apply(counts,1,function(x) sum(x<10))<ncol(counts)*0.9,]
nrow(counts4wgcna)
ncol(counts4wgcna)
write.csv(counts4wgcna, file="../../outputs/transcriptomics/ssid/deseq2/counts4wgcna.csv")

# importing a design .csv file
design = read.csv("../../data/transcriptomics/ssid/design_ssid.csv", head=TRUE)
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
outs=c(34,76)
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ site+treat)
dds$site <- factor(dds$site, levels = c("Emerald", "Rainbow", "Star", "MacN"))
dds$treat <- factor(dds$treat, levels = c("CC", "LC", "CH", "LH"))

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="../../outputs/transcriptomics/ssid/deseq2/initial.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,2],design[,11],sep=".")
# renames the column names
colnames(vsd)=snames

save(vsd,design,file="../../outputs/transcriptomics/ssid/deseq2/vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ site+treat)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=TRUE)
# vsd.wg=assay(rlog(wg), blind=TRUE)
head(vsd.wg)
colnames(vsd.wg)=snames
save(vsd.wg,design,file="../../outputs/transcriptomics/ssid/deseq2/data4wgcna.RData")


#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("../../outputs/transcriptomics/ssid/deseq2/vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="../../outputs/transcriptomics/ssid/deseq2/heatmap_ssid_host.pdf", width=15, height=15)
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
write.csv(dds.pcoa$values, file = "../../outputs/transcriptomics/ssid/deseq2/PCoA_variance.csv")

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues 
pdf(file="../../outputs/transcriptomics/ssid/deseq2/PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's

# plotting PCoA by site and treatment
pdf(file="../../outputs/transcriptomics/ssid/deseq2/PCoA_ssid_host.pdf", width=12, height=6)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))],pch=c(6,2,25,17)[as.numeric(as.factor(conditions$treat))], bg= c("#018571","#80cdc1","#dfc27d","#a6611a")[as.numeric(as.factor(conditions$site))], xlab="Coordinate 1 (16.6%)", ylab="Coordinate 2 (4.6%)", main="Site")
ordiellipse(scores, conditions$site, label=F, col=c("#018571","#80cdc1","#dfc27d","#a6611a"))
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), fill = c("#018571","#80cdc1","#dfc27d","#a6611a"), bty="n")
legend("topleft", legend=c("CC", "CH", "LC", "LH"), pch=c(6,2,25,17), bty="n", pt.bg = c(NA,NA,"black",NA))
plot(scores[,1], scores[,2],col=c("#92c5de","#f4a582","#0571b0","#ca0020")[as.numeric(as.factor(conditions$treat))],pch=c(15,0,1,16)[as.numeric((as.factor(conditions$site)))], xlab="Coordinate 1 (16.6%)", ylab="Coordinate 2 (4.6%)", main="Treatment")
ordiellipse(scores, conditions$treat, label=F, col=c("#92c5de","#f4a582","#0571b0","#ca0020"))
legend("topleft", legend=c("CC", "CH", "LC", "LH"), fill = c("#92c5de","#f4a582","#0571b0","#ca0020"), bty="n")
legend("topright", legend=c("Emerald", "Rainbow", "Star", "MacN"), pch=c(15,0,1,16), bty="n")
dev.off()

# neighbor-joining tree of samples (based on significant PCo's):
pdf(file="../../outputs/transcriptomics/ssid/deseq2/PCoA_tree.pdf", width=10, height=10)
tre=nj(dist(scores[,1:4]))
plot(tre,cex=0.8)
dev.off()

# formal analysis of variance in distance matricies: 
set.seed(23485)
ad=adonis2(t(vsd)~site*treat, data=conditions, method="manhattan", by = "terms", permutations=1e6)
ad
write.csv(ad, file = "../../outputs/transcriptomics/ssid/deseq2/PERMANOVA_output.csv")

# creating pie chart to represent ANOVA results
cols=c("blue","orange","grey80")
pdf(file="../../outputs/transcriptomics/ssid/deseq2/PERMANOVA_pie.pdf", width=6, height=6)
pie(ad$R2[1:3],labels=row.names(ad)[1:4],col=cols,main="site vs treatment")
dev.off()


#### DESEQ ####

# with multi-factor, multi-level design - using LRT
load("../../outputs/transcriptomics/ssid/deseq2/initial.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
dds=DESeq(dds, parallel=TRUE)

# model for the effect of treatment: (>2 factor levels => LRT)
dds$treat <- factor(dds$treat, levels = c("CC", "LC", "CH", "LH"))
dds_treat=DESeq(dds,test="LRT",reduced=~site, parallel=TRUE)

# saving all models
save(dds,dds_treat,file="../../outputs/transcriptomics/ssid/deseq2/realModels.RData")


#### DEGs and CONTRASTS ####

load("../../outputs/transcriptomics/ssid/deseq2/realModels.RData")
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
     degs_LC_CC, degs_CH_CC, degs_LH_CC, degs_CH_LC, degs_LH_CH, degs_LH_LC, file="../../outputs/transcriptomics/ssid/deseq2/pvals.RData")


#### VENN DIAGRAMS ####

load("../../outputs/transcriptomics/ssid/deseq2/pvals.RData")
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
pdf(file="../../outputs/transcriptomics/ssid/deseq2/Venn_site.pdf", height=10, width=12)
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
pdf(file="../../outputs/transcriptomics/ssid/deseq2/Venn_treatment.pdf", height=12, width=12)
grid.draw(venn_treat)
dev.off()


#### GO/KOG EXPORT ####

load("../../outputs/transcriptomics/ssid/deseq2/realModels.RData")
load("../../outputs/transcriptomics/ssid/deseq2/pvals.RData")

# site contrasts

# signed log p-values: -log(padj)* direction:
source=Rainbow_Emerald[!is.na(Rainbow_Emerald$pvalue),]
Rainbow_Emerald.p=data.frame("gene"=row.names(source))
Rainbow_Emerald.p$lpv=-log(source[,"padj"],10)
Rainbow_Emerald.p$lpv[source$stat<0]=Rainbow_Emerald.p$lpv[source$stat<0]*-1
head(Rainbow_Emerald.p)
write.csv(Rainbow_Emerald.p,file="../../outputs/transcriptomics/ssid/deseq2/Rainbow_Emerald_lpv.csv",row.names=F,quote=F)
save(Rainbow_Emerald.p,file="../../outputs/transcriptomics/ssid/deseq2/Rainbow_Emerald_lpv.RData")

source=Star_Emerald[!is.na(Star_Emerald$pvalue),]
Star_Emerald.p=data.frame("gene"=row.names(source))
Star_Emerald.p$lpv=-log(source[,"padj"],10)
Star_Emerald.p$lpv[source$stat<0]=Star_Emerald.p$lpv[source$stat<0]*-1
head(Star_Emerald.p)
write.csv(Star_Emerald.p,file="../../outputs/transcriptomics/ssid/deseq2/Star_Emerald_lpv.csv",row.names=F,quote=F)
save(Star_Emerald.p,file="../../outputs/transcriptomics/ssid/deseq2/Star_Emerald_lpv.RData")

source=MacN_Emerald[!is.na(MacN_Emerald$pvalue),]
MacN_Emerald.p=data.frame("gene"=row.names(source))
MacN_Emerald.p$lpv=-log(source[,"padj"],10)
MacN_Emerald.p$lpv[source$stat<0]=MacN_Emerald.p$lpv[source$stat<0]*-1
head(MacN_Emerald.p)
write.csv(MacN_Emerald.p,file="../../outputs/transcriptomics/ssid/deseq2/MacN_Emerald_lpv.csv",row.names=F,quote=F)
save(MacN_Emerald.p,file="../../outputs/transcriptomics/ssid/deseq2/MacN_Emerald_lpv.RData")

source=Star_Rainbow[!is.na(Star_Rainbow$pvalue),]
Star_Rainbow.p=data.frame("gene"=row.names(source))
Star_Rainbow.p$lpv=-log(source[,"padj"],10)
Star_Rainbow.p$lpv[source$stat<0]=Star_Rainbow.p$lpv[source$stat<0]*-1
head(Star_Rainbow.p)
write.csv(Star_Rainbow.p,file="../../outputs/transcriptomics/ssid/deseq2/Star_Rainbow_lpv.csv",row.names=F,quote=F)
save(Star_Rainbow.p,file="../../outputs/transcriptomics/ssid/deseq2/Star_Rainbow_lpv.RData")

source=MacN_Rainbow[!is.na(MacN_Rainbow$pvalue),]
MacN_Rainbow.p=data.frame("gene"=row.names(source))
MacN_Rainbow.p$lpv=-log(source[,"padj"],10)
MacN_Rainbow.p$lpv[source$stat<0]=MacN_Rainbow.p$lpv[source$stat<0]*-1
head(MacN_Rainbow.p)
write.csv(MacN_Rainbow.p,file="../../outputs/transcriptomics/ssid/deseq2/MacN_Rainbow_lpv.csv",row.names=F,quote=F)
save(MacN_Rainbow.p,file="../../outputs/transcriptomics/ssid/deseq2/MacN_Rainbow_lpv.RData")

source=MacN_Star[!is.na(MacN_Star$pvalue),]
MacN_Star.p=data.frame("gene"=row.names(source))
MacN_Star.p$lpv=-log(source[,"padj"],10)
MacN_Star.p$lpv[source$stat<0]=MacN_Star.p$lpv[source$stat<0]*-1
head(MacN_Star.p)
write.csv(MacN_Star.p,file="../../outputs/transcriptomics/ssid/deseq2/MacN_Star_lpv.csv",row.names=F,quote=F)
save(MacN_Star.p,file="../../outputs/transcriptomics/ssid/deseq2/MacN_Star_lpv.RData")

# treatment contrasts

# signed log p-values: -log(padj)* direction:
source=LC_CC[!is.na(LC_CC$pvalue),]
LC_CC.p=data.frame("gene"=row.names(source))
LC_CC.p$lpv=-log(source[,"padj"],10)
LC_CC.p$lpv[source$stat<0]=LC_CC.p$lpv[source$stat<0]*-1
head(LC_CC.p)
write.csv(LC_CC.p,file="../../outputs/transcriptomics/ssid/deseq2/LC_CC_lpv.csv",row.names=F,quote=F)
save(LC_CC.p,file="../../outputs/transcriptomics/ssid/deseq2/LC_CC_lpv.RData")

source=CH_CC[!is.na(CH_CC$pvalue),]
CH_CC.p=data.frame("gene"=row.names(source))
CH_CC.p$lpv=-log(source[,"padj"],10)
CH_CC.p$lpv[source$stat<0]=CH_CC.p$lpv[source$stat<0]*-1
head(CH_CC.p)
write.csv(CH_CC.p,file="../../outputs/transcriptomics/ssid/deseq2/CH_CC_lpv.csv",row.names=F,quote=F)
save(CH_CC.p,file="../../outputs/transcriptomics/ssid/deseq2/CH_CC_lpv.RData")

source=LH_CC[!is.na(LH_CC$pvalue),]
LH_CC.p=data.frame("gene"=row.names(source))
LH_CC.p$lpv=-log(source[,"padj"],10)
LH_CC.p$lpv[source$stat<0]=LH_CC.p$lpv[source$stat<0]*-1
head(LH_CC.p)
write.csv(LH_CC.p,file="../../outputs/transcriptomics/ssid/deseq2/LH_CC_lpv.csv",row.names=F,quote=F)
save(LH_CC.p,file="../../outputs/transcriptomics/ssid/deseq2/LH_CC_lpv.RData")

source=CH_LC[!is.na(CH_LC$pvalue),]
CH_LC.p=data.frame("gene"=row.names(source))
CH_LC.p$lpv=-log(source[,"padj"],10)
CH_LC.p$lpv[source$stat<0]=CH_LC.p$lpv[source$stat<0]*-1
head(CH_LC.p)
write.csv(CH_LC.p,file="../../outputs/transcriptomics/ssid/deseq2/CH_LC_lpv.csv",row.names=F,quote=F)
save(CH_LC.p,file="../../outputs/transcriptomics/ssid/deseq2/CH_LC_lpv.RData")

source=LH_LC[!is.na(LH_LC$pvalue),]
LH_LC.p=data.frame("gene"=row.names(source))
LH_LC.p$lpv=-log(source[,"padj"],10)
LH_LC.p$lpv[source$stat<0]=LH_LC.p$lpv[source$stat<0]*-1
head(LH_LC.p)
write.csv(LH_LC.p,file="../../outputs/transcriptomics/ssid/deseq2/LH_LC_lpv.csv",row.names=F,quote=F)
save(LH_LC.p,file="../../outputs/transcriptomics/ssid/deseq2/LH_LC_lpv.RData")

source=LH_CH[!is.na(LH_CH$pvalue),]
LH_CH.p=data.frame("gene"=row.names(source))
LH_CH.p$lpv=-log(source[,"padj"],10)
LH_CH.p$lpv[source$stat<0]=LH_CH.p$lpv[source$stat<0]*-1
head(LH_CH.p)
write.csv(LH_CH.p,file="../../outputs/transcriptomics/ssid/deseq2/LH_CH_lpv.csv",row.names=F,quote=F)
save(LH_CH.p,file="../../outputs/transcriptomics/ssid/deseq2/LH_CH_lpv.RData")

save(Rainbow_Emerald.p, Star_Emerald.p, MacN_Emerald.p, Star_Rainbow.p, MacN_Rainbow.p, MacN_Star.p, 
     LC_CC.p, CH_CC.p, LH_CC.p, CH_LC.p, LH_CH.p, LH_LC.p, file="../../outputs/transcriptomics/ssid/deseq2/exports.RData")


#### DEG MATCHING ####

library(DESeq2)
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
load("../../outputs/transcriptomics/ssid/deseq2/exports.RData")

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
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="stress_control", .before="gene") -> commongenes_treatment

# exporting all DEGs matching across stress vs control treatments
write.csv(commongenes_treatment, file="../../outputs/transcriptomics/ssid/deseq2/commongenes_treatment.csv")

# pairwise treatment comparisons
LC_CC.p %>%
  dplyr::inner_join(CH_CC.p, by = "gene", suffix = c(".LC_CC", ".CH_CC")) %>%
  filter(abs(lpv.LC_CC) >= 1 & abs(lpv.CH_CC) >= 1) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> bleaching_oa

LC_CC.p %>%
  dplyr::inner_join(LH_CC.p, by = "gene", suffix = c(".LC_CC", ".LH_CC")) %>%
  filter(abs(lpv.LC_CC) >= 1 & abs(lpv.LH_CC) >= 1) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> oableaching_oa

CH_CC.p %>%
  dplyr::inner_join(LH_CC.p, by = "gene", suffix = c(".CH_CC", ".LH_CC")) %>%
  filter(abs(lpv.CH_CC) >= 1 & abs(lpv.LH_CC) >= 1) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> bleaching_oableaching

# joining all matching DEGs into a single dataframe
commongenes_pairwise_treatment <- bind_rows(bleaching_oa,oableaching_oa,bleaching_oableaching)
write.csv(commongenes_pairwise_treatment, file="../../outputs/transcriptomics/ssid/deseq2/commongenes_pairwise_treatment.csv")

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
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="urban_reef", .before="gene") -> commongenes_site

# exporting all DEGs matching across stress vs control treatments
write.csv(commongenes_site, file="../../outputs/transcriptomics/ssid/deseq2/commongenes_site.csv")

# pairwise site comparisons
Star_Emerald.p %>%
  dplyr::inner_join(MacN_Emerald.p, by = "gene", suffix =c(".Star_Emerald", ".MacN_Emerald")) %>%
  filter(abs(lpv.Star_Emerald) >= 1 & abs(lpv.MacN_Emerald) >= 1) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> urban_Emerald

Star_Rainbow.p %>%
  dplyr::inner_join(MacN_Rainbow.p, by = "gene", suffix =c(".Star_Rainbow", ".MacN_Rainbow")) %>%
  filter(abs(lpv.Star_Rainbow) >= 1 & abs(lpv.MacN_Rainbow) >= 1) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Siderastrea-siderea-annotated-transcriptome/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> urban_Rainbow


# joining all matching DEGs into a single dataframe
commongenes_pairwise_site <- bind_rows(urban_Emerald,urban_Rainbow)
write.csv(commongenes_pairwise_site, file="../../outputs/transcriptomics/ssid/deseq2/commongenes_pairwise_site.csv")

#### CORRELATION PLOT HELPERS ####

# --- Unicode-safe symbols (avoid warnings on some devices) ---
use_unicode <- l10n_info()[["UTF-8"]] && capabilities("cairo")
LTE   <- if (use_unicode) "≤" else "<="
MINUS <- if (use_unicode) "−" else "-"

# --- Parameters ---
inf_cap  <- 20
REL_NAME <- "Relationship"
SIG_NAME <- "p Value"

# Discordant filter knobs (used by both dataframe sections)
delta_cut_dd <- 2.0
ratio_cut_dd <- 100
high_cut_dd  <- 2.0   # ~ p <= 0.01
low_cut_dd   <- 1.3   # ~ p <= 0.05
min_sig      <- 1

# --- Relationship colors (colorblind-friendly) ---
REL_LEVELS <- c("Direct", "Inverse")
rel_cols   <- setNames(c("#009E73", "#D55E00"), REL_LEVELS)

# --- Significance tiers (shared by all plots) ---
SIG_LEVELS <- c("p > 0.1", "p ≤ 0.05", "p ≤ 1e−6")
sig_alpha  <- setNames(c(0.15, 0.90, 0.90), SIG_LEVELS)
sig_size   <- setNames(c(1.0,  1.2,  1.2),  SIG_LEVELS)
sig_stroke <- setNames(c(0.00, 0.20, 0.00), SIG_LEVELS)

# --- Shared helper functions ---
# Normalize any two named lpv columns into lpv_x / lpv_y, cap infinities, drop non-finite rows
clean_lpv_df <- function(df, col_x, col_y, cap = inf_cap) {
  if (is.null(df) || !nrow(df)) stop("Input data frame is empty.", call. = FALSE)
  if (!all(c(col_x, col_y) %in% names(df)))
    stop(sprintf("Columns not found: %s",
                 paste(setdiff(c(col_x, col_y), names(df)), collapse = ", ")), call. = FALSE)
  d <- df
  d$lpv_x <- as.numeric(d[[col_x]])
  d$lpv_y <- as.numeric(d[[col_y]])
  d$lpv_x[ is.infinite(d$lpv_x) & d$lpv_x > 0 ] <- cap
  d$lpv_x[ is.infinite(d$lpv_x) & d$lpv_x < 0 ] <- -cap
  d$lpv_y[ is.infinite(d$lpv_y) & d$lpv_y > 0 ] <- cap
  d$lpv_y[ is.infinite(d$lpv_y) & d$lpv_y < 0 ] <- -cap
  d <- d[is.finite(d$lpv_x) & is.finite(d$lpv_y), , drop = FALSE]
  d
}

safe_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- sum(complete.cases(x, y))
  if (n >= 3 && sd(x) > 0 && sd(y) > 0) {
    ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
    r  <- unname(ct$estimate); p <- ct$p.value
    list(r = r, p = p, n = n)
  } else list(r = NA_real_, p = NA_real_, n = n)
}

fmt_num <- function(x) {
  if (is.na(x)) "NA" else if (abs(x) < 1e-3) format(x, digits = 2, scientific = TRUE) else sprintf("%.3f", x)
}

sig_category3 <- function(lpv_x, lpv_y, levels_vec = SIG_LEVELS) {
  lv <- pmin(abs(lpv_x), abs(lpv_y))
  p  <- 10^(-lv)
  out <- dplyr::case_when(
    p <= 1e-6 ~ levels_vec[3],
    p <= 0.05 ~ levels_vec[2],
    TRUE      ~ levels_vec[1]
  )
  factor(out, levels = levels_vec)
}

compute_global_limits <- function(comparisons, min_span = 0.6, pad_frac = 0.05) {
  vals <- unlist(lapply(comparisons, function(entry) {
    d <- clean_lpv_df(entry$df, entry$col_x, entry$col_y)
    d <- d[pmax(abs(d$lpv_x), abs(d$lpv_y)) >= 1.0, ]
    c(d$lpv_x, d$lpv_y)
  }))
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs)) max_abs <- min_span / 2
  lim <- max(max_abs * (1 + pad_frac), min_span / 2)
  c(-lim, lim)
}

make_panel <- function(df, title, col_x, col_y, x_lab, y_lab, global_lim) {
  dat <- clean_lpv_df(df, col_x, col_y)
  dat <- subset(dat, pmax(abs(dat$lpv_x), abs(dat$lpv_y)) >= 1.0)

  same <- (dat$lpv_x >= 0 & dat$lpv_y >= 0) | (dat$lpv_x <= 0 & dat$lpv_y <= 0)
  dat$rel_class <- factor(ifelse(same, "Direct", "Inverse"), levels = REL_LEVELS)
  dat$sig_cat   <- sig_category3(dat$lpv_x, dat$lpv_y, levels_vec = SIG_LEVELS)

  x_limits <- global_lim
  y_limits <- global_lim

  span       <- diff(global_lim)
  jitter_amt <- span * 0.006
  set.seed(42)
  dat$x_jit <- pmin(pmax(jitter(dat$lpv_x, amount = jitter_amt), global_lim[1]), global_lim[2])
  dat$y_jit <- pmin(pmax(jitter(dat$lpv_y, amount = jitter_amt), global_lim[1]), global_lim[2])

  st <- safe_stats(dat$lpv_x, dat$lpv_y)
  lab_text <- if (is.na(st$r)) {
    sprintf("n = %d", st$n)
  } else {
    sprintf("r = %.2f, p = %s, n = %d", st$r, fmt_num(st$p), st$n)
  }

  ref_levels <- c("Direct 1:1", "Inverse 1:1")

  ggplot(dat, aes(x = lpv_x, y = lpv_y)) +
    geom_hline(yintercept = 0, color = "grey85") +
    geom_vline(xintercept = 0, color = "grey85") +

    geom_point(
      aes(x = x_jit, y = y_jit, color = rel_class, fill = rel_class,
          alpha = sig_cat, stroke = sig_cat, size = sig_cat),
      shape = 21
    ) +

    geom_point(
      data = subset(dat, sig_cat == SIG_LEVELS[3]),
      aes(x = x_jit, y = y_jit, fill = rel_class, alpha = sig_cat, size = sig_cat),
      color = "grey15", shape = 21, stroke = 0.40
    ) +

    geom_point(
      data = data.frame(sig_cat = factor(SIG_LEVELS, levels = SIG_LEVELS)),
      mapping = aes(x = 0, y = 0, alpha = sig_cat),
      inherit.aes = FALSE, shape = 21, size = 0, stroke = 0, fill = NA, color = NA,
      show.legend = TRUE, na.rm = TRUE
    ) +

    geom_abline(intercept = 0, slope =  1, color = "firebrick", linewidth = 1,
                linetype = "dashed", alpha = 0.5, show.legend = FALSE) +
    geom_abline(intercept = 0, slope = -1, color = "firebrick", linewidth = 1,
                linetype = "dotted", alpha = 0.5, show.legend = FALSE) +

    geom_segment(
      data = data.frame(x = 0, y = 0, xend = 0, yend = 0,
                        ref = factor(ref_levels, levels = ref_levels)),
      mapping = aes(x = x, y = y, xend = xend, yend = yend, linetype = ref),
      inherit.aes = FALSE, color = "firebrick", linewidth = 1, alpha = 0.5,
      show.legend = TRUE
    ) +

    annotate("text", x = x_limits[1], y = y_limits[2], label = lab_text,
             hjust = 0, vjust = 1.1, size = 4.2) +

    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +

    scale_color_manual(values = rel_cols, breaks = REL_LEVELS, limits = REL_LEVELS,
                       drop = FALSE, guide = "none") +
    scale_discrete_manual("size", values = sig_size, breaks = SIG_LEVELS, limits = SIG_LEVELS,
                          guide = "none") +
    scale_fill_manual(values = rel_cols, breaks = REL_LEVELS, limits = REL_LEVELS,
                      drop = FALSE, name = REL_NAME) +
    scale_alpha_manual(values = sig_alpha, breaks = SIG_LEVELS, limits = SIG_LEVELS,
                       drop = FALSE, name = SIG_NAME) +
    scale_discrete_manual("stroke", values = sig_stroke,
                          breaks = SIG_LEVELS, limits = SIG_LEVELS, guide = "none") +
    scale_linetype_manual(
      name   = "Reference",
      values = c("Direct 1:1" = "dashed", "Inverse 1:1" = "dotted")
    ) +
    guides(
      fill  = guide_legend(order = 1,
                           override.aes = list(shape = 21, size = 1.6, linetype = 0,
                                               stroke = c(0.45, 0.45))),
      alpha = guide_legend(order = 2,
                           override.aes = list(shape = 21, fill = "grey50",
                                               color = c("grey30", "grey30", "grey15"),
                                               size = 1.6, linetype = 0,
                                               stroke = c(0.00, 0.20, 0.40))),
      linetype = guide_legend(order = 3,
                              override.aes = list(color = "firebrick", linewidth = 1, alpha = 0.5))
    ) +
    labs(
      x = str_wrap(paste0(x_lab, " (signed ", MINUS, "log10 p)"), width = 28),
      y = str_wrap(paste0(y_lab, " (signed ", MINUS, "log10 p)"), width = 28),
      title = str_wrap(title, width = 32)
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "right",
      axis.line = element_blank(),
      plot.margin = margin(8, 20, 8, 8)
    )
}

classify_for_filters <- function(df, col_x, col_y) {
  d <- clean_lpv_df(df, col_x, col_y)
  same <- (d$lpv_x >= 0 & d$lpv_y >= 0) | (d$lpv_x <= 0 & d$lpv_y <= 0)
  d$rel_class    <- ifelse(same, "Direct", "Inverse")
  d$abs_x        <- abs(d$lpv_x)
  d$abs_y        <- abs(d$lpv_y)
  d$delta_log10p <- abs(d$abs_x - d$abs_y)
  d$p_x          <- 10^(-d$abs_x)
  d$p_y          <- 10^(-d$abs_y)
  d$p_ratio      <- pmax(d$p_x, d$p_y) / pmin(d$p_x, d$p_y)
  d
}

rbind_align <- function(...) {
  dfs  <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(x) {
    if (is.null(x) || !nrow(x)) {
      return(as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols)))
    }
    miss <- setdiff(cols, names(x))
    for (m in miss) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, dfs2)
}

build_panels <- function(comparisons, ncol = 3) {
  global_lim <- compute_global_limits(comparisons)
  panels <- mapply(
    FUN      = function(title, entry) make_panel(entry$df, title, entry$col_x, entry$col_y,
                                                 entry$x_lab, entry$y_lab, global_lim),
    title    = names(comparisons),
    entry    = comparisons,
    SIMPLIFY = FALSE
  )
  nrow_grid <- ceiling(length(panels) / ncol)
  for (i in seq_along(panels)) {
    row <- ceiling(i / ncol)
    col <- i - (row - 1) * ncol
    if (col != 1) panels[[i]] <- panels[[i]] + theme(axis.title.y = element_blank())
    if (row != nrow_grid) panels[[i]] <- panels[[i]] + theme(axis.title.x = element_blank())
  }
  one_legend  <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
  legend_plot <- cowplot::ggdraw(one_legend)
  panels_noleg <- lapply(panels, function(p) p + theme(legend.position = "none"))
  grid  <- wrap_plots(panels_noleg, ncol = ncol)
  panel <- grid | legend_plot
  panel + plot_layout(widths = c(1, 0.30))
}

filter_discordant <- function(comparisons, group_col) {
  filter_one <- function(entry, group_name) {
    d         <- classify_for_filters(entry$df, entry$col_x, entry$col_y)
    meets_min <- pmax(d$abs_x, d$abs_y) >= min_sig

    inverse <- d[d$rel_class == "Inverse" & meets_min, , drop = FALSE]
    if (nrow(inverse)) {
      inverse$filter_type           <- "Inverse"
      inverse$relationship_category <- "Inverse"
      inverse$discord_direction     <- NA_character_
      inverse$delta_cut_used        <- NA_real_
      inverse$high_cut_used         <- NA_real_
      inverse$low_cut_used          <- NA_real_
      inverse$ratio_cut_used        <- NA_real_
      inverse[[group_col]]          <- group_name
    }

    strong_x_weak_y <- d$abs_x >= high_cut_dd & d$abs_y <= low_cut_dd
    strong_y_weak_x <- d$abs_y >= high_cut_dd & d$abs_x <= low_cut_dd
    cross_tier <- strong_x_weak_y | strong_y_weak_x
    dd_idx <- d$rel_class == "Direct" &
      d$delta_log10p >= delta_cut_dd &
      d$p_ratio      >= ratio_cut_dd &
      cross_tier & meets_min

    direct_discordant <- d[dd_idx, , drop = FALSE]
    if (nrow(direct_discordant)) {
      so <- direct_discordant$abs_x >= high_cut_dd & direct_discordant$abs_y <= low_cut_dd
      direct_discordant$filter_type           <- "DirectDiscordant"
      direct_discordant$relationship_category <- "Direct (discordant significance)"
      direct_discordant$discord_direction     <- ifelse(so,
                                                        paste0(entry$x_lab, " >> ", entry$y_lab),
                                                        paste0(entry$y_lab, " >> ", entry$x_lab))
      direct_discordant$delta_cut_used        <- delta_cut_dd
      direct_discordant$high_cut_used         <- high_cut_dd
      direct_discordant$low_cut_used          <- low_cut_dd
      direct_discordant$ratio_cut_used        <- ratio_cut_dd
      direct_discordant[[group_col]]          <- group_name
    }

    rbind_align(inverse, direct_discordant)
  }

  do.call(rbind_align, lapply(names(comparisons), function(nm) filter_one(comparisons[[nm]], nm)))
}


#### CORRELATION PLOTS TREATMENT ####

treatment_comparisons <- list(
  "Bleaching vs OA" = list(
    df    = bleaching_oa,
    col_x = "lpv.LC_CC", col_y = "lpv.CH_CC",
    x_lab = "OA",         y_lab = "Bleaching"),
  "OA + Bleaching vs OA" = list(
    df    = oableaching_oa,
    col_x = "lpv.LC_CC", col_y = "lpv.LH_CC",
    x_lab = "OA",         y_lab = "OA + Bleaching"),
  "OA + Bleaching vs Bleaching" = list(
    df    = bleaching_oableaching,
    col_x = "lpv.CH_CC", col_y = "lpv.LH_CC",
    x_lab = "Bleaching",  y_lab = "OA + Bleaching")
)

panel_treatment <- build_panels(treatment_comparisons, ncol = 3)
print(panel_treatment)
ggsave("../../outputs/transcriptomics/ssid/deseq2/commongenes correlation treatment.pdf",
       panel_treatment, width = 14, height = 5.5)


#### CORRELATION DATAFRAME TREATMENT ####

treatment_filtered_genes <- filter_discordant(treatment_comparisons, group_col = "treatment")

if (nrow(treatment_filtered_genes) > 0) {
  cat("Rows total:", nrow(treatment_filtered_genes), "\n")
  print(table(treatment_filtered_genes$relationship_category, useNA = "ifany"))
  print(table(treatment_filtered_genes$treatment, useNA = "ifany"))
} else {
  warning("treatment_filtered_genes is empty with current strict thresholds.")
}
dplyr::count(treatment_filtered_genes, treatment, relationship_category)
write.csv(treatment_filtered_genes,
          file = "../../outputs/transcriptomics/ssid/deseq2/commongenes correlation discordant treatment.csv")


#### CORRELATION PLOTS SITE ####

site_comparisons <- list(
  "Urban Sites vs Emerald Reef" = list(
    df    = urban_Emerald,
    col_x = "lpv.Star_Emerald",  col_y = "lpv.MacN_Emerald",
    x_lab = "Star Island vs Emerald Reef",
    y_lab = "MacArthur North vs Emerald Reef"),
  "Urban Sites vs Rainbow Reef" = list(
    df    = urban_Rainbow,
    col_x = "lpv.Star_Rainbow",  col_y = "lpv.MacN_Rainbow",
    x_lab = "Star Island vs Rainbow Reef",
    y_lab = "MacArthur North vs Rainbow Reef")
)

panel_site <- build_panels(site_comparisons, ncol = 2)
print(panel_site)
ggsave("../../outputs/transcriptomics/ssid/deseq2/commongenes correlation site.pdf",
       panel_site, width = 9, height = 5.5)


#### CORRELATION DATAFRAME SITE ####

site_filtered_genes <- filter_discordant(site_comparisons, group_col = "site")

if (nrow(site_filtered_genes) > 0) {
  cat("Rows total:", nrow(site_filtered_genes), "\n")
  print(table(site_filtered_genes$relationship_category, useNA = "ifany"))
  print(table(site_filtered_genes$site, useNA = "ifany"))
} else {
  warning("site_filtered_genes is empty with current strict thresholds.")
}
dplyr::count(site_filtered_genes, site, relationship_category)
write.csv(site_filtered_genes,
          file = "../../outputs/transcriptomics/ssid/deseq2/commongenes correlation discordant site.csv")


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
source("../../scripts/transcriptomics/uniHeatmap.R")

# Treatment
# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(commongenes_treatment_heatmap$gene, commongenes_treatment_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1
pdf(file="../../outputs/transcriptomics/ssid/deseq2/heatmap_treatment_p0.1.pdf", height=25, width=50)
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
gene_names2 <- as.data.frame(cbind(commongenes_site_heatmap$gene, commongenes_site_heatmap$annot))

# p < 0.1
pdf(file="../../outputs/transcriptomics/ssid/deseq2/heatmap_site_p0.1.pdf", height=24, width=48)
uniHeatmap(vsd=vsd_site,gene.names=gene_names2,
           metric=-(abs(commongenes_site_heatmap$lpv.MacN_Emerald)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
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
ggsave("../../outputs/transcriptomics/ssid/deseq2/KOG_treatment.pdf", plot= KOG_sum, width=10, height=6, units="in", dpi=300)

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
ggsave("../../outputs/transcriptomics/ssid/deseq2/KOG_site.pdf", plot= KOG_sum, width=14, height=6, units="in", dpi=300)

