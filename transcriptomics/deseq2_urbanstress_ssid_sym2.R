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
counts = read.table("allcounts_sym2.txt")

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
design = read.csv("design_ssid.csv", head=TRUE)
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

# CANNOT PROCEED - not enough genes