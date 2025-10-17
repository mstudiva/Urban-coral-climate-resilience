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
library(tidyverse)
library(stringr)

pair_panels <- list(
  "LC_CC vs CH_CC" =
    make_pairpanel_cols(commongenes_treatment, "lpv.LC_CC", "lpv.CH_CC") +
    labs(
      x = str_wrap("Acidified + Ambient vs Contemporary + Ambient", width = 40),
      y = str_wrap("Contemporary + Bleaching vs Contemporary + Ambient", width = 40),
      title = NULL
    ) +
    theme(plot.title = element_blank()),

  "LC_CC vs LH_CC" =
    make_pairpanel_cols(commongenes_treatment, "lpv.LC_CC", "lpv.LH_CC") +
    labs(
      x = str_wrap("Acidified + Ambient vs Contemporary + Ambient", width = 40),
      y = str_wrap("Acidified + Bleaching vs Contemporary + Ambient", width = 40),
      title = NULL
    ) +
    theme(plot.title = element_blank()),

  "CH_CC vs LH_CC" =
    make_pairpanel_cols(commongenes_treatment, "lpv.CH_CC", "lpv.LH_CC") +
    labs(
      x = str_wrap("Contemporary + Bleaching vs Contemporary + Ambient", width = 40),
      y = str_wrap("Acidified + Bleaching vs Contemporary + Ambient", width = 40),
      title = NULL
    ) +
    theme(plot.title = element_blank())
)

# Hide redundant Y-axis title for 2nd and 3rd panels
for (i in seq_along(pair_panels)) {
  if (i != 1) pair_panels[[i]] <- pair_panels[[i]] + theme(axis.title.y = element_blank())
}

# Single legend & layout
one_legend  <- cowplot::get_legend(pair_panels[[1]] + theme(legend.position = "right"))
legend_plot <- cowplot::ggdraw(one_legend)
panels_noleg <- lapply(pair_panels, function(p) p + theme(legend.position = "none"))

grid  <- patchwork::wrap_plots(panels_noleg, ncol = 3)
panel <- grid | legend_plot
panel <- panel + patchwork::plot_layout(widths = c(1, 0.30))

print(panel)
ggsave("pairwise_contrast_panels_labels_wrapped.pdf", panel, width = 14, height = 5.5)


#read in counts
counts = read.table("../../../../raw/ofav/allcounts_host.txt")

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
design = read.csv("../../../../raw/design_ofav.csv", head=TRUE)
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

source=LH_LC[!is.na(LH_LC$pvalue),]
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
library(patchwork)
load("exports.RData")

# This section of code does several things: 1) join -log10(pval) across treatment comparisons, 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3) adds gene annotations, and 4) then pulls on corresponding KOG classes

# Treatment
# stress treatments (LC, CH, and LH) versus control treatment (CC)
LC_CC.p %>%
  dplyr::inner_join(CH_CC.p, by = "gene", suffix = c(".LC_CC", ".CH_CC")) %>%
  dplyr::inner_join(LH_CC.p %>% dplyr::rename(lpv.LH_CC = lpv), by = "gene") %>%
  filter(abs(lpv.LC_CC) >= 1 & abs(lpv.CH_CC) >= 1 & abs(lpv.LH_CC) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="stress_control", .before="gene") -> commongenes_treatment

# exporting all DEGs matching across stress vs control treatments
write.csv(commongenes_treatment, file="commongenes_treatment.csv")

# pairwise treatment comparisons
LC_CC.p %>%
  dplyr::inner_join(CH_CC.p, by = "gene", suffix = c(".LC_CC", ".CH_CC")) %>%
  filter(abs(lpv.LC_CC) >= 1 & abs(lpv.CH_CC) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> bleaching_oa

LC_CC.p %>%
  dplyr::inner_join(LH_CC.p, by = "gene", suffix = c(".LC_CC", ".LH_CC")) %>%
  filter(abs(lpv.LC_CC) >= 1 & abs(lpv.LH_CC) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> oableaching_oa

CH_CC.p %>%
  dplyr::inner_join(LH_CC.p, by = "gene", suffix = c(".CH_CC", ".LH_CC")) %>%
  filter(abs(lpv.CH_CC) >= 1 & abs(lpv.LH_CC) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> bleaching_oableaching

# joining all matching DEGs into a single dataframe
commongenes_pairwise_treatment <- bind_rows(bleaching_oa,oableaching_oa,bleaching_oableaching)
write.csv(commongenes_pairwise_treatment, file="commongenes_pairwise_treatment.csv")

# Site
# This section of code does several things: 1) join -log10(pval) across site comparisons, 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3) adds gene annotations, and 4) then pulls on corresponding KOG classes

# urban sites (Star and MacN) versus reef sites (Emerald and Rainbow)
Star_Emerald.p %>%
  dplyr::inner_join(MacN_Emerald.p, by = "gene", suffix =c(".Star_Emerald", ".MacN_Emerald")) %>%
  dplyr::inner_join(Star_Rainbow.p %>% dplyr::rename(lpv.Star_Rainbow = lpv), by = "gene") %>%
  dplyr::inner_join(MacN_Rainbow.p %>% dplyr::rename(lpv.MacN_Rainbow = lpv), by = "gene") %>%
  filter(abs(lpv.Star_Emerald) >= 1 & abs(lpv.MacN_Emerald) >= 1 & abs(lpv.Star_Rainbow) >= 1 & abs(lpv.MacN_Rainbow) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="urban_reef", .before="gene") -> commongenes_site

# exporting all DEGs matching across stress vs control treatments
write.csv(commongenes_site, file="commongenes_site.csv")

# pairwise site comparisons
Star_Emerald.p %>%
  dplyr::inner_join(MacN_Emerald.p, by = "gene", suffix =c(".Star_Emerald", ".MacN_Emerald")) %>%
  filter(abs(lpv.Star_Emerald) >= 1 & abs(lpv.MacN_Emerald) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> urban_Emerald

Star_Rainbow.p %>%
  dplyr::inner_join(MacN_Rainbow.p, by = "gene", suffix =c(".Star_Rainbow", ".MacN_Rainbow")) %>%
  filter(abs(lpv.Star_Rainbow) >= 1 & abs(lpv.MacN_Rainbow) >= 1) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) -> urban_Rainbow

# joining all matching DEGs into a single dataframe
commongenes_pairwise_site <- bind_rows(urban_Emerald,urban_Rainbow)
write.csv(commongenes_pairwise_site, file="commongenes_pairwise_site.csv")

#### CORRELATION PLOTS TREATMENT ####
# --- Unicode-safe symbols (avoid warnings on some devices) ---
use_unicode <- l10n_info()[["UTF-8"]] && capabilities("cairo")
LTE   <- if (use_unicode) "\u2264" else "<="
MINUS <- if (use_unicode) "\u2212" else "-"

# --- Parameters (tweak as needed) ---
sig_cut1_global <- 1.3  # ~ p = 0.05
sig_cut2_global <- 2.0  # ~ p = 0.01
sig_cut3_global <- 6.0  # ~ p = 1e-6  (raise to 4.0 for ~1e-4)
inf_cap          <- 20   # cap |−log10 p| when p == 0 (Inf)
REL_NAME <- "Relationship"
SIG_NAME <- "p Value"

# --- Relationship colors (colorblind-friendly) ---
REL_LEVELS <- c("Direct", "Inverse")
rel_cols   <- setNames(c("#009E73", "#D55E00"), REL_LEVELS)  # green, orange

# --- Helpers (generic; no species assumptions) ---
clean_lpv_df <- function(df, col_x, col_y, cap = inf_cap) {
  if (is.null(df) || !nrow(df)) stop("Input data frame is empty.", call. = FALSE)
  if (!all(c(col_x, col_y) %in% names(df))) {
    stop(sprintf("Columns not found: %s",
                 paste(setdiff(c(col_x, col_y), names(df)), collapse = ", ")), call. = FALSE)
  }
  d <- df
  d$lpv_x <- as.numeric(d[[col_x]])
  d$lpv_y <- as.numeric(d[[col_y]])
  
  # Cap infinities from p=0
  d$lpv_x[ is.infinite(d$lpv_x) & d$lpv_x > 0 ] <- cap
  d$lpv_x[ is.infinite(d$lpv_x) & d$lpv_x < 0 ] <- -cap
  d$lpv_y[ is.infinite(d$lpv_y) & d$lpv_y > 0 ] <- cap
  d$lpv_y[ is.infinite(d$lpv_y) & d$lpv_y < 0 ] <- -cap
  
  # Drop non-finite rows
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

# --- p-value label helpers (legend shows ACTUAL p ranges) ---
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) sprintf("%.0e", p) else sprintf("%.3f", p)
}

pvalue_labels <- function(c1, c2, c3) {
  v <- sort(c(c1, c2, c3)); c1 <- v[1]; c2 <- v[2]; c3 <- v[3]
  p1 <- 10^(-c1); p2 <- 10^(-c2); p3 <- 10^(-c3)
  c(
    sprintf("p > %s",           fmt_p(p1)),
    sprintf("%s < p %s %s",     fmt_p(p2), LTE, fmt_p(p1)),
    sprintf("%s < p %s %s",     fmt_p(p3), LTE, fmt_p(p2)),
    sprintf("p %s %s",          LTE, fmt_p(p3))
  )
}

# Significance tiers (alpha, light -> dark)
SIG_LEVELS <- pvalue_labels(sig_cut1_global, sig_cut2_global, sig_cut3_global)
sig_alpha  <- setNames(c(0.15, 0.35, 0.65, 1.00), SIG_LEVELS)

sig_category4 <- function(lpv_x, lpv_y, c1, c2, c3, levels_vec = SIG_LEVELS) {
  v <- sort(c(c1, c2, c3)); c1 <- v[1]; c2 <- v[2]; c3 <- v[3]
  lv <- pmin(abs(lpv_x), abs(lpv_y))  # min significance across the pair in |−log10 p|
  out <- cut(
    lv,
    breaks = c(-Inf, c1, c2, c3, Inf),
    labels = pvalue_labels(c1, c2, c3),
    include.lowest = TRUE, right = TRUE
  )
  factor(out, levels = levels_vec)
}

# --- One-panel constructor (explicit columns & labels; dashed/dotted refs) ---
make_panel <- function(df, title, col_x, col_y, x_lab, y_lab,
                       c1 = sig_cut1_global, c2 = sig_cut2_global, c3 = sig_cut3_global) {
  dat <- clean_lpv_df(df, col_x, col_y)
  
  # Relationship (direct/inverse) based on signs
  same <- (dat$lpv_x >= 0 & dat$lpv_y >= 0) | (dat$lpv_x <= 0 & dat$lpv_y <= 0)
  dat$rel_class <- factor(ifelse(same, "Direct", "Inverse"), levels = REL_LEVELS)
  
  # Significance tiers
  dat$sig_cat <- sig_category4(dat$lpv_x, dat$lpv_y, c1, c2, c3, levels_vec = SIG_LEVELS)
  
  # --- Data-driven square limits (not forced around 0), with min span ---
  x_min <- min(dat$lpv_x, na.rm = TRUE); x_max <- max(dat$lpv_x, na.rm = TRUE)
  y_min <- min(dat$lpv_y, na.rm = TRUE); y_max <- max(dat$lpv_y, na.rm = TRUE)
  x_center <- (x_min + x_max) / 2
  y_center <- (y_min + y_max) / 2
  
  # base spans (handle single-point cases)
  x_span0 <- max(x_max - x_min, 1e-8)
  y_span0 <- max(y_max - y_min, 1e-8)
  
  # jitter + padding
  jitter_w <- x_span0 * 0.003
  jitter_h <- y_span0 * 0.003
  pad_x <- max(0.04 * x_span0, 6 * jitter_w)
  pad_y <- max(0.04 * y_span0, 6 * jitter_h)
  
  # padded spans
  x_span <- x_span0 + 2 * pad_x
  y_span <- y_span0 + 2 * pad_y
  
  # enforce a minimum square span so tiny panels don’t collapse
  min_span <- 0.6
  final_span <- max(x_span, y_span, min_span)
  
  x_limits <- c(x_center - final_span / 2, x_center + final_span / 2)
  y_limits <- c(y_center - final_span / 2, y_center + final_span / 2)
  
  # Stats label
  st <- safe_stats(dat$lpv_x, dat$lpv_y)
  fmt_num <- function(x) if (is.na(x)) "NA" else if (abs(x) < 1e-3) format(x, digits = 2, scientific = TRUE) else sprintf("%.3f", x)
  lab_text <- if (is.na(st$r)) sprintf("n = %d", st$n) else sprintf("r = %.2f, p = %s, n = %d", st$r, fmt_num(st$p), st$n)
  
  # Legend labels for line types
  ref_levels <- c("Direct 1:1", "Inverse 1:1")
  
  ggplot(dat, aes(x = lpv_x, y = lpv_y)) +
    # Draw origin axes only if 0 is inside the panel
    (if (y_limits[1] < 0 && y_limits[2] > 0) geom_hline(yintercept = 0, color = "grey85") else NULL) +
    (if (x_limits[1] < 0 && x_limits[2] > 0) geom_vline(xintercept = 0, color = "grey85") else NULL) +
    
    # Points (linetype=NA so legend circles stay clean)
    geom_point(
      aes(color = rel_class, fill = rel_class, alpha = sig_cat),
      shape = 21, size = 1.6, stroke = 0.45, linetype = NA,
      position = position_jitter(width = final_span * 0.003, height = final_span * 0.003)
    ) +
    
    # Alpha legend trainer (invisible)
    geom_point(
      data = data.frame(sig_cat = factor(SIG_LEVELS, levels = SIG_LEVELS)),
      mapping = aes(x = 0, y = 0, alpha = sig_cat),
      inherit.aes = FALSE, shape = 21, size = 0, stroke = 0, fill = NA, color = NA,
      show.legend = TRUE
    ) +
    
    # Reference lines (drawn; legend fed by trainer below)
    geom_abline(intercept = 0, slope =  1, color = "firebrick", linewidth = 1,
                linetype = "dashed", alpha = 0.5, show.legend = FALSE) +
    geom_abline(intercept = 0, slope = -1, color = "firebrick", linewidth = 1,
                linetype = "dotted", alpha = 0.5, show.legend = FALSE) +
    
    # Line-type legend trainer (zero-length segments create legend keys)
    geom_segment(
      data = data.frame(x = 0, y = 0, xend = 0, yend = 0,
                        ref = factor(ref_levels, levels = ref_levels)),
      mapping = aes(x = x, y = y, xend = xend, yend = yend, linetype = ref),
      inherit.aes = FALSE, color = "firebrick", linewidth = 1, alpha = 0.5,
      show.legend = TRUE
    ) +
    
    # Stats label (top-left inside the frame)
    annotate("text", x = x_limits[1], y = y_limits[2], label = lab_text,
             hjust = 0, vjust = 1.1, size = 4.2) +
    
    # Clip inside panel with square limits
    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +
    
    # Legends & scales
    scale_color_manual(values = rel_cols, breaks = REL_LEVELS, limits = REL_LEVELS, drop = FALSE, guide = "none") +
    scale_fill_manual(values = rel_cols,  breaks = REL_LEVELS, limits = REL_LEVELS, drop = FALSE, name = REL_NAME) +
    scale_alpha_manual(values = sig_alpha, breaks = SIG_LEVELS, limits = SIG_LEVELS, drop = FALSE, name = SIG_NAME) +
    scale_linetype_manual(
      name   = "Reference",
      values = c("Direct 1:1" = "dashed", "Inverse 1:1" = "dotted")
    ) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 1.6, stroke = 0.45, linetype = 0)),
      alpha = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey50", color = "grey30", size = 1.6, stroke = 0.45, linetype = 0)),
      linetype = guide_legend(order = 3, override.aes = list(color = "firebrick", linewidth = 1, alpha = 0.5))
    ) +
    labs(
      x = paste0(x_lab, " (signed ", MINUS, "log10 p)"),
      y = paste0(y_lab, " (signed ", MINUS, "log10 p)"),
      title = NULL
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

# --- Build 3 panels (explicit columns & labels) ---
# Data frames and exact columns:
# bleaching_oa        : x = lpv.LC_CC, y = lpv.CH_CC
# oableaching_oa      : x = lpv.LC_CC, y = lpv.LH_CC
# bleaching_oableaching: x = lpv.CH_CC, y = lpv.LH_CC

comparisons <- list(
  list(title = "Bleaching vs OA",
       df    = bleaching_oa,
       col_x = "lpv.LC_CC", col_y = "lpv.CH_CC",
       x_lab = "OA",            y_lab = "Bleaching"),
  list(title = "OA+Bleaching vs OA",
       df    = oableaching_oa,
       col_x = "lpv.LC_CC", col_y = "lpv.LH_CC",
       x_lab = "OA",            y_lab = "OA + Bleaching"),
  list(title = "Bleaching vs OA+Bleaching",
       df    = bleaching_oableaching,
       col_x = "lpv.CH_CC", col_y = "lpv.LH_CC",
       x_lab = "Bleaching",     y_lab = "OA + Bleaching")
)

panels <- lapply(
  comparisons,
  function(cmp) make_panel(
    df     = cmp$df,
    title  = cmp$title,
    col_x  = cmp$col_x,
    col_y  = cmp$col_y,
    x_lab  = cmp$x_lab,
    y_lab  = cmp$y_lab,
    c1 = sig_cut1_global, c2 = sig_cut2_global, c3 = sig_cut3_global
  )
)

# --- Hide redundant axis TITLES only (keep ticks & numbers) for a 1x3 layout ---
ncol_grid <- 3
nrow_grid <- ceiling(length(panels) / ncol_grid)

for (i in seq_along(panels)) {
  row <- ceiling(i / ncol_grid)
  # Keep ALL y-axis titles: do NOT blank axis.title.y
  if (row != nrow_grid) {
    panels[[i]] <- panels[[i]] + theme(axis.title.x = element_blank())
  }
}

# --- Single legend workflow ---
one_legend  <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
legend_plot <- cowplot::ggdraw(one_legend)

panels_noleg <- lapply(panels, function(p) p + theme(legend.position = "none"))

grid  <- wrap_plots(panels_noleg, ncol = ncol_grid)
panel <- grid | legend_plot
panel <- panel + plot_layout(widths = c(1, 0.30))

# --- Show & Save ---
print(panel)
ggsave("commongenes correlation treatment.pdf", panel, width = 14, height = 5.5)
# ggsave("commongenes correlation treatment.png", panel, width = 14, height = 5.5, dpi = 600)


#### CORRELATION DATAFRAME TREATMENT ####
# --- knobs ---
delta_cut_dd <- 2.0   # |Δ −log10 p| >= 2  => ≥100x p-value difference
ratio_cut_dd <- 100   # p-ratio >= 100
high_cut_dd  <- sig_cut2_global   # strong (e.g., 2.0)
low_cut_dd   <- sig_cut1_global   # weak   (e.g., 1.3)
min_sig      <- 0.0

# --- comparison definitions (columns & pretty axis names) ---
comparisons_df <- list(
  list(name = "Bleaching vs OA",
       df   = bleaching_oa,
       col_x = "lpv.LC_CC", col_y = "lpv.CH_CC",
       x_nm = "OA",         y_nm = "Bleaching"),
  list(name = "OA+Bleaching vs OA",
       df   = oableaching_oa,
       col_x = "lpv.LC_CC", col_y = "lpv.LH_CC",
       x_nm = "OA",         y_nm = "OA + Bleaching"),
  list(name = "Bleaching vs OA+Bleaching",
       df   = bleaching_oableaching,
       col_x = "lpv.CH_CC", col_y = "lpv.LH_CC",
       x_nm = "Bleaching",  y_nm = "OA + Bleaching")
)

# --- align & rbind helper (handles empties) ---
rbind_align <- function(...) {
  dfs  <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(x) {
    if (is.null(x) || !nrow(x)) {
      as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
    } else {
      miss <- setdiff(cols, names(x))
      for (m in miss) x[[m]] <- NA
      x[, cols, drop = FALSE]
    }
  })
  do.call(rbind, dfs2)
}

# --- derived fields for filtering (uses your clean_lpv_df(df, col_x, col_y)) ---
classify_for_filters_xy <- function(df, col_x, col_y, comp_name, x_nm, y_nm) {
  d <- clean_lpv_df(df, col_x, col_y)  # creates lpv_x, lpv_y
  
  same <- (d$lpv_x >= 0 & d$lpv_y >= 0) | (d$lpv_x <= 0 & d$lpv_y <= 0)
  d$rel_class    <- ifelse(same, "Direct", "Inverse")
  
  d$abs_x        <- abs(d$lpv_x)
  d$abs_y        <- abs(d$lpv_y)
  d$delta_log10p <- abs(d$abs_x - d$abs_y)
  d$p_x          <- 10^(-d$abs_x)
  d$p_y          <- 10^(-d$abs_y)
  d$p_ratio      <- pmax(d$p_x, d$p_y) / pmin(d$p_x, d$p_y)
  
  d$comparison   <- comp_name
  d$x_axis_name  <- x_nm
  d$y_axis_name  <- y_nm
  d$x_col        <- col_x
  d$y_col        <- col_y
  d
}

# --- filter one comparison into inverse + strict direct-discordant ---
filter_one_xy <- function(df, comp_name, col_x, col_y, x_nm, y_nm) {
  d <- classify_for_filters_xy(df, col_x, col_y, comp_name, x_nm, y_nm)
  meets_min <- (pmax(d$abs_x, d$abs_y) >= min_sig)
  
  # 1) ALL inverse
  inverse <- d[d$rel_class == "Inverse" & meets_min, , drop = FALSE]
  if (nrow(inverse)) {
    inverse$filter_type           <- "inverse"
    inverse$relationship_category <- "Inverse"
    inverse$discord_direction     <- NA_character_
    inverse$delta_cut_used        <- NA_real_
    inverse$high_cut_used         <- NA_real_
    inverse$low_cut_used          <- NA_real_
    inverse$ratio_cut_used        <- NA_real_
  }
  
  # 2) STRICT direct-discordant
  strong_x_weak_y_full <- d$abs_x >= high_cut_dd & d$abs_y <= low_cut_dd
  strong_y_weak_x_full <- d$abs_y >= high_cut_dd & d$abs_x <= low_cut_dd
  cross_tier_full      <- strong_x_weak_y_full | strong_y_weak_x_full
  
  dd_idx <- d$rel_class == "Direct" &
    d$delta_log10p >= delta_cut_dd &
    d$p_ratio      >= ratio_cut_dd &
    cross_tier_full &
    meets_min
  
  discordant <- d[dd_idx, , drop = FALSE]
  if (nrow(discordant)) {
    # compute direction INSIDE the subset to match lengths
    strong_x_weak_y_sub <- discordant$abs_x >= high_cut_dd & discordant$abs_y <= low_cut_dd
    discordant$filter_type           <- "discordant"
    discordant$relationship_category <- "Direct (discordant significance)"
    discordant$discord_direction     <- ifelse(strong_x_weak_y_sub,
                                               paste0(x_nm, " >> ", y_nm),
                                               paste0(y_nm, " >> ", x_nm))
    discordant$delta_cut_used        <- delta_cut_dd
    discordant$high_cut_used         <- high_cut_dd
    discordant$low_cut_used          <- low_cut_dd
    discordant$ratio_cut_used        <- ratio_cut_dd
  }
  
  list(inverse = inverse, discordant = discordant)
}

# --- run filters across all comparisons ---
res_list <- lapply(
  comparisons_df,
  function(cmp) filter_one_xy(
    df       = cmp$df,
    comp_name= cmp$name,
    col_x    = cmp$col_x, col_y = cmp$col_y,
    x_nm     = cmp$x_nm,  y_nm  = cmp$y_nm
  )
)

# --- build ONE MASTER DATAFRAME with filter_type column ---
inverse_list    <- lapply(res_list, `[[`, "inverse")
discordant_list <- lapply(res_list, `[[`, "discordant")

inverse_discordant_master <- do.call(
  rbind_align,
  c(inverse_list, discordant_list)
)

# --- Build row labels and write ONE CSV (ID first, comparison second) ---
if (exists("inverse_discordant_master") && nrow(inverse_discordant_master)) {
  df <- inverse_discordant_master
  
  # sequential index within each comparison
  idx_within_comp <- ave(seq_len(nrow(df)), df$comparison, FUN = seq_along)
  
  # row label: "<comparison>_<n>"
  df$row_label <- paste(df$comparison, idx_within_comp, sep = "_")
  
  # put row_label first, comparison second; keep the rest as-is
  lead2 <- c("row_label", "comparison")
  rest  <- setdiff(names(df), lead2)
  df_out <- df[, c(lead2, rest), drop = FALSE]
  
  # sanity prints (optional)
  cat("Total rows in master:", nrow(df_out), "\n")
  print(table(df_out$filter_type, useNA = "ifany"))
  print(table(df_out$comparison,  useNA = "ifany"))
  
  # write ONE CSV (this is the only write)
  write.csv(df_out, file = "commongenes correlation discordant treatment.csv", row.names = FALSE)
} else {
  warning("inverse_discordant_master is missing or empty; no CSV written.")
}


#### CORRELATION PLOTS SITE ####
# --- Unicode-safe symbols (avoid warnings on some devices) ---
use_unicode <- l10n_info()[["UTF-8"]] && capabilities("cairo")
LTE   <- if (use_unicode) "\u2264" else "<="
MINUS <- if (use_unicode) "\u2212" else "-"

# --- Parameters (tweak as needed) ---
sig_cut1_global <- 1.3  # ~ p = 0.05
sig_cut2_global <- 2.0  # ~ p = 0.01
sig_cut3_global <- 6.0  # ~ p = 1e-6  (raise to 4.0 for ~1e-4)
inf_cap          <- 20   # cap |−log10 p| when p == 0 (Inf)
REL_NAME <- "Relationship"
SIG_NAME <- "p Value"

# --- Relationship colors (colorblind-friendly) ---
REL_LEVELS <- c("Direct", "Inverse")
rel_cols   <- setNames(c("#009E73", "#D55E00"), REL_LEVELS)  # green, orange

# --- Helpers (generic; no species assumptions) ---
clean_lpv_df <- function(df, col_x, col_y, cap = inf_cap) {
  if (is.null(df) || !nrow(df)) stop("Input data frame is empty.", call. = FALSE)
  if (!all(c(col_x, col_y) %in% names(df))) {
    stop(sprintf("Columns not found: %s",
                 paste(setdiff(c(col_x, col_y), names(df)), collapse = ", ")), call. = FALSE)
  }
  d <- df
  d$lpv_x <- as.numeric(d[[col_x]])
  d$lpv_y <- as.numeric(d[[col_y]])
  
  # Cap infinities from p=0
  d$lpv_x[ is.infinite(d$lpv_x) & d$lpv_x > 0 ] <- cap
  d$lpv_x[ is.infinite(d$lpv_x) & d$lpv_x < 0 ] <- -cap
  d$lpv_y[ is.infinite(d$lpv_y) & d$lpv_y > 0 ] <- cap
  d$lpv_y[ is.infinite(d$lpv_y) & d$lpv_y < 0 ] <- -cap
  
  # Drop non-finite rows
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

# --- p-value label helpers (legend shows ACTUAL p ranges) ---
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) sprintf("%.0e", p) else sprintf("%.3f", p)
}

pvalue_labels <- function(c1, c2, c3) {
  v <- sort(c(c1, c2, c3)); c1 <- v[1]; c2 <- v[2]; c3 <- v[3]
  p1 <- 10^(-c1); p2 <- 10^(-c2); p3 <- 10^(-c3)
  c(
    sprintf("p > %s",           fmt_p(p1)),
    sprintf("%s < p %s %s",     fmt_p(p2), LTE, fmt_p(p1)),
    sprintf("%s < p %s %s",     fmt_p(p3), LTE, fmt_p(p2)),
    sprintf("p %s %s",          LTE, fmt_p(p3))
  )
}

# Significance tiers (alpha, light -> dark)
SIG_LEVELS <- pvalue_labels(sig_cut1_global, sig_cut2_global, sig_cut3_global)
sig_alpha  <- setNames(c(0.15, 0.35, 0.65, 1.00), SIG_LEVELS)

sig_category4 <- function(lpv_x, lpv_y, c1, c2, c3, levels_vec = SIG_LEVELS) {
  v <- sort(c(c1, c2, c3)); c1 <- v[1]; c2 <- v[2]; c3 <- v[3]
  lv <- pmin(abs(lpv_x), abs(lpv_y))  # min significance across the pair in |−log10 p|
  out <- cut(
    lv,
    breaks = c(-Inf, c1, c2, c3, Inf),
    labels = pvalue_labels(c1, c2, c3),
    include.lowest = TRUE, right = TRUE
  )
  factor(out, levels = levels_vec)
}

# --- One-panel constructor (explicit columns & labels; dashed/dotted refs) ---
make_panel <- function(df, title, col_x, col_y, x_lab, y_lab,
                       c1 = sig_cut1_global, c2 = sig_cut2_global, c3 = sig_cut3_global) {
  dat <- clean_lpv_df(df, col_x, col_y)
  
  # Relationship (direct/inverse) based on signs
  same <- (dat$lpv_x >= 0 & dat$lpv_y >= 0) | (dat$lpv_x <= 0 & dat$lpv_y <= 0)
  dat$rel_class <- factor(ifelse(same, "Direct", "Inverse"), levels = REL_LEVELS)
  
  # Significance tiers
  dat$sig_cat <- sig_category4(dat$lpv_x, dat$lpv_y, c1, c2, c3, levels_vec = SIG_LEVELS)
  
  # --- Data-driven square limits (not forced around 0), with min span ---
  x_min <- min(dat$lpv_x, na.rm = TRUE); x_max <- max(dat$lpv_x, na.rm = TRUE)
  y_min <- min(dat$lpv_y, na.rm = TRUE); y_max <- max(dat$lpv_y, na.rm = TRUE)
  x_center <- (x_min + x_max) / 2
  y_center <- (y_min + y_max) / 2
  
  # base spans (handle single-point cases)
  x_span0 <- max(x_max - x_min, 1e-8)
  y_span0 <- max(y_max - y_min, 1e-8)
  
  # jitter + padding
  jitter_w <- x_span0 * 0.003
  jitter_h <- y_span0 * 0.003
  pad_x <- max(0.04 * x_span0, 6 * jitter_w)
  pad_y <- max(0.04 * y_span0, 6 * jitter_h)
  
  # padded spans
  x_span <- x_span0 + 2 * pad_x
  y_span <- y_span0 + 2 * pad_y
  
  # enforce a minimum square span so tiny panels don’t collapse
  min_span <- 0.6
  final_span <- max(x_span, y_span, min_span)
  
  x_limits <- c(x_center - final_span / 2, x_center + final_span / 2)
  y_limits <- c(y_center - final_span / 2, y_center + final_span / 2)
  
  # Stats label
  st <- safe_stats(dat$lpv_x, dat$lpv_y)
  fmt_num <- function(x) if (is.na(x)) "NA" else if (abs(x) < 1e-3) format(x, digits = 2, scientific = TRUE) else sprintf("%.3f", x)
  lab_text <- if (is.na(st$r)) sprintf("n = %d", st$n) else sprintf("r = %.2f, p = %s, n = %d", st$r, fmt_num(st$p), st$n)
  
  # Legend labels for line types
  ref_levels <- c("Direct 1:1", "Inverse 1:1")
  
  ggplot(dat, aes(x = lpv_x, y = lpv_y)) +
    # Draw origin axes only if 0 is inside the panel
    (if (y_limits[1] < 0 && y_limits[2] > 0) geom_hline(yintercept = 0, color = "grey85") else NULL) +
    (if (x_limits[1] < 0 && x_limits[2] > 0) geom_vline(xintercept = 0, color = "grey85") else NULL) +
    
    # Points (linetype=NA so legend circles stay clean)
    geom_point(
      aes(color = rel_class, fill = rel_class, alpha = sig_cat),
      shape = 21, size = 1.6, stroke = 0.45, linetype = NA,
      position = position_jitter(width = final_span * 0.003, height = final_span * 0.003)
    ) +
    
    # Alpha legend trainer (invisible)
    geom_point(
      data = data.frame(sig_cat = factor(SIG_LEVELS, levels = SIG_LEVELS)),
      mapping = aes(x = 0, y = 0, alpha = sig_cat),
      inherit.aes = FALSE, shape = 21, size = 0, stroke = 0, fill = NA, color = NA,
      show.legend = TRUE
    ) +
    
    # Reference lines (drawn; legend fed by trainer below)
    geom_abline(intercept = 0, slope =  1, color = "firebrick", linewidth = 1,
                linetype = "dashed", alpha = 0.5, show.legend = FALSE) +
    geom_abline(intercept = 0, slope = -1, color = "firebrick", linewidth = 1,
                linetype = "dotted", alpha = 0.5, show.legend = FALSE) +
    
    # Line-type legend trainer (zero-length segments create legend keys)
    geom_segment(
      data = data.frame(x = 0, y = 0, xend = 0, yend = 0,
                        ref = factor(ref_levels, levels = ref_levels)),
      mapping = aes(x = x, y = y, xend = xend, yend = yend, linetype = ref),
      inherit.aes = FALSE, color = "firebrick", linewidth = 1, alpha = 0.5,
      show.legend = TRUE
    ) +
    
    # Stats label (top-left inside the frame)
    annotate("text", x = x_limits[1], y = y_limits[2], label = lab_text,
             hjust = 0, vjust = 1.1, size = 4.2) +
    
    # Clip inside panel with square limits
    coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE, clip = "on") +
    
    # Legends & scales
    scale_color_manual(values = rel_cols, breaks = REL_LEVELS, limits = REL_LEVELS, drop = FALSE, guide = "none") +
    scale_fill_manual(values = rel_cols,  breaks = REL_LEVELS, limits = REL_LEVELS, drop = FALSE, name = REL_NAME) +
    scale_alpha_manual(values = sig_alpha, breaks = SIG_LEVELS, limits = SIG_LEVELS, drop = FALSE, name = SIG_NAME) +
    scale_linetype_manual(
      name   = "Reference",
      values = c("Direct 1:1" = "dashed", "Inverse 1:1" = "dotted")
    ) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 1.6, stroke = 0.45, linetype = 0)),
      alpha = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey50", color = "grey30", size = 1.6, stroke = 0.45, linetype = 0)),
      linetype = guide_legend(order = 3, override.aes = list(color = "firebrick", linewidth = 1, alpha = 0.5))
    ) +
    labs(
      x = str_wrap(paste0(x_lab, " (signed ", MINUS, "log10 p)"), width = 28),
      y = str_wrap(paste0(y_lab, " (signed ", MINUS, "log10 p)"), width = 28),
      title = NULL
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

# --- Build 2 panels (explicit columns & labels) ---
# Star_reef:   x = lpv.Star_Emerald,  y = lpv.Star_Rainbow
# MacN_reef:   x = lpv.MacN_Emerald,  y = lpv.MacN_Rainbow

comparisons <- list(
  list(title = "Star and MacN vs Emerald",
       df    = urban_Emerald,
       col_x = "lpv.Star_Emerald", col_y = "lpv.MacN_Emerald",
       x_lab = "Star Island vs Emerald Reef",
       y_lab = "MacArthur North vs Emerald Reef"),
  list(title = "Star and MacN vs Rainbow",
       df    = urban_Rainbow,
       col_x = "lpv.Star_Rainbow", col_y = "lpv.MacN_Rainbow",
       x_lab = "Star Island vs Rainbow Reef",
       y_lab = "MacArthur North vs Rainbow Reef")
)

panels <- lapply(
  comparisons,
  function(cmp) make_panel(
    df     = cmp$df,
    title  = cmp$title,
    col_x  = cmp$col_x,
    col_y  = cmp$col_y,
    x_lab  = cmp$x_lab,
    y_lab  = cmp$y_lab,
    c1 = sig_cut1_global, c2 = sig_cut2_global, c3 = sig_cut3_global
  )
)

# --- Keep ALL y-axis titles; hide x titles only on non-bottom rows (1 row => keep all) ---
ncol_grid <- 2
nrow_grid <- ceiling(length(panels) / ncol_grid)
for (i in seq_along(panels)) {
  row <- ceiling(i / ncol_grid)
  if (row != nrow_grid) {
    panels[[i]] <- panels[[i]] + theme(axis.title.x = element_blank())
  }
}

# --- Single legend workflow ---
one_legend  <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
legend_plot <- cowplot::ggdraw(one_legend)

panels_noleg <- lapply(panels, function(p) p + theme(legend.position = "none"))

grid  <- wrap_plots(panels_noleg, ncol = ncol_grid)
panel <- grid | legend_plot
panel <- panel + plot_layout(widths = c(1, 0.30))

# --- Show & Save ---
print(panel)
ggsave("commongenes correlation site.pdf", panel, width = 12, height = 5.5)
# ggsave("commongenes correlation site.png", panel, width = 12, height = 5.5, dpi = 600)


#### CORRELATION DATAFRAME SITE ####
# --- knobs ---
delta_cut_dd <- 2.0   # |Δ −log10 p| >= 2  => ≥100x p-value difference
ratio_cut_dd <- 100   # p-ratio >= 100
high_cut_dd  <- sig_cut2_global   # strong (e.g., 2.0)
low_cut_dd   <- sig_cut1_global   # weak   (e.g., 1.3)
min_sig      <- 0.0

# --- site comparison definitions (columns & pretty axis names) ---
comparisons_df <- list(
  list(name = "Urban Sites vs Emerald Reef",
       df   = urban_Emerald,
       col_x = "lpv.Star_Emerald",  col_y = "lpv.MacN_Emerald",
       x_nm = "Star Island vs Emerald Reef",
       y_nm = "MacArthur North vs Emerald Reef"),
  list(name = "Urban Sites vs Rainbow Reef",
       df   = urban_Rainbow,
       col_x = "lpv.Star_Rainbow",  col_y = "lpv.MacN_Rainbow",
       x_nm = "Star Island vs Rainbow Reef",
       y_nm = "MacArthur North vs Rainbow Reef")
)

# --- align & rbind helper (handles empties) ---
rbind_align <- function(...) {
  dfs  <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(x) {
    if (is.null(x) || !nrow(x)) {
      as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
    } else {
      miss <- setdiff(cols, names(x))
      for (m in miss) x[[m]] <- NA
      x[, cols, drop = FALSE]
    }
  })
  do.call(rbind, dfs2)
}

# --- derived fields for filtering (uses clean_lpv_df(df, col_x, col_y)) ---
classify_for_filters_xy <- function(df, col_x, col_y, comp_name, x_nm, y_nm) {
  d <- clean_lpv_df(df, col_x, col_y)  # creates lpv_x, lpv_y
  
  # relationship (by sign)
  same <- (d$lpv_x >= 0 & d$lpv_y >= 0) | (d$lpv_x <= 0 & d$lpv_y <= 0)
  d$rel_class    <- ifelse(same, "Direct", "Inverse")
  
  # magnitudes and p-values
  d$abs_x        <- abs(d$lpv_x)
  d$abs_y        <- abs(d$lpv_y)
  d$delta_log10p <- abs(d$abs_x - d$abs_y)
  d$p_x          <- 10^(-d$abs_x)
  d$p_y          <- 10^(-d$abs_y)
  d$p_ratio      <- pmax(d$p_x, d$p_y) / pmin(d$p_x, d$p_y)
  
  # meta
  d$comparison   <- comp_name
  d$x_axis_name  <- x_nm
  d$y_axis_name  <- y_nm
  d$x_col        <- col_x
  d$y_col        <- col_y
  d
}

# --- filter one comparison into inverse + strict direct-discordant ---
filter_one_xy <- function(df, comp_name, col_x, col_y, x_nm, y_nm) {
  d <- classify_for_filters_xy(df, col_x, col_y, comp_name, x_nm, y_nm)
  meets_min <- (pmax(d$abs_x, d$abs_y) >= min_sig)
  
  # 1) ALL inverse
  inverse <- d[d$rel_class == "Inverse" & meets_min, , drop = FALSE]
  if (nrow(inverse)) {
    inverse$filter_type           <- "inverse"
    inverse$relationship_category <- "Inverse"
    inverse$discord_direction     <- NA_character_
    inverse$delta_cut_used        <- NA_real_
    inverse$high_cut_used         <- NA_real_
    inverse$low_cut_used          <- NA_real_
    inverse$ratio_cut_used        <- NA_real_
  }
  
  # 2) STRICT direct-discordant
  strong_x_weak_y_full <- d$abs_x >= high_cut_dd & d$abs_y <= low_cut_dd
  strong_y_weak_x_full <- d$abs_y >= high_cut_dd & d$abs_x <= low_cut_dd
  cross_tier_full      <- strong_x_weak_y_full | strong_y_weak_x_full
  
  dd_idx <- d$rel_class == "Direct" &
    d$delta_log10p >= delta_cut_dd &
    d$p_ratio      >= ratio_cut_dd &
    cross_tier_full &
    meets_min
  
  discordant <- d[dd_idx, , drop = FALSE]
  if (nrow(discordant)) {
    # compute direction INSIDE the subset to match lengths
    strong_x_weak_y_sub <- discordant$abs_x >= high_cut_dd & discordant$abs_y <= low_cut_dd
    discordant$filter_type           <- "discordant"
    discordant$relationship_category <- "Direct (discordant significance)"
    discordant$discord_direction     <- ifelse(strong_x_weak_y_sub,
                                               paste0(x_nm, " >> ", y_nm),
                                               paste0(y_nm, " >> ", x_nm))
    discordant$delta_cut_used        <- delta_cut_dd
    discordant$high_cut_used         <- high_cut_dd
    discordant$low_cut_used          <- low_cut_dd
    discordant$ratio_cut_used        <- ratio_cut_dd
  }
  
  list(inverse = inverse, discordant = discordant)
}

# --- run filters across both site comparisons ---
res_list <- lapply(
  comparisons_df,
  function(cmp) filter_one_xy(
    df       = cmp$df,
    comp_name= cmp$name,
    col_x    = cmp$col_x, col_y = cmp$col_y,
    x_nm     = cmp$x_nm,  y_nm  = cmp$y_nm
  )
)

# --- build ONE MASTER DATAFRAME with filter_type column ---
inverse_list    <- lapply(res_list, `[[`, "inverse")
discordant_list <- lapply(res_list, `[[`, "discordant")

inverse_discordant_master <- do.call(
  rbind_align,
  c(inverse_list, discordant_list)
)

# --- Build row labels and write ONE CSV (ID first, comparison second) ---
if (exists("inverse_discordant_master") && nrow(inverse_discordant_master)) {
  df <- inverse_discordant_master
  
  # sequential index within each comparison
  idx_within_comp <- ave(seq_len(nrow(df)), df$comparison, FUN = seq_along)
  
  # row label: "<comparison>_<n>"
  df$row_label <- paste(df$comparison, idx_within_comp, sep = "_")
  
  # put row_label first, comparison second; keep the rest as-is
  lead2 <- c("row_label", "comparison")
  rest  <- setdiff(names(df), lead2)
  df_out <- df[, c(lead2, rest), drop = FALSE]
  
  # sanity prints (optional)
  cat("Total rows in master:", nrow(df_out), "\n")
  print(table(df_out$filter_type, useNA = "ifany"))
  print(table(df_out$comparison,  useNA = "ifany"))
  
  # write ONE CSV
  write.csv(df_out, file = "commongenes correlation discordant site.csv", row.names = FALSE)
} else {
  warning("inverse_discordant_master is missing or empty; no CSV written.")
}


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