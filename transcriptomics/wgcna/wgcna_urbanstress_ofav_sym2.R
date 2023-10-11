#### PACKAGES ####

# installing WGCNA:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.17")
# BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("flashClust")
# install.packages("WGCNA",dependencies=TRUE)
# repos="http://cran.us.r-project.org"
# run these above commands once, then comment out

# always run these before running any of the following script chunks
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
# ?WGCNA


#### DATA IMPORT and TRAITS ####

# importing data generated from DESeq2 script
lnames=load("../../../DESeq2/ofav/sym2/data4wgcna.RData")
lnames # "vsd.wg"  "design" # log-transformed variance-stabilized gene expression, and table or experimental conditions
datt=t(vsd.wg)
ncol(datt)
nrow(datt)

head(design)
str(design)

# assembling table of traits

# coding experimental factors as binary (0/1, yes/no)
Emerald=as.numeric(design$site=="Emerald") 
Rainbow=as.numeric(design$site=="Rainbow") 
Star=as.numeric(design$site=="Star") 
MacN=as.numeric(design$site=="MacN") 
lowpH=as.numeric(design$ph=="lowpH") 
hightemp=as.numeric(design$temp=="hightemp") 

# combining with numerical traits (physiological data)
traits <- cbind(Emerald, Rainbow, Star, MacN, lowpH, hightemp, design[c(13:21)])
traits


#### OUTLIER DETECTION ####

# identifies outlier genes
gsg = goodSamplesGenes(datt, verbose = 3);
gsg$allOK #if TRUE, no outlier genes

# calculates mean expression per array, then the number of missing values per array
meanExpressionByArray=apply( datt,1,mean, na.rm=T)
NumberMissingByArray=apply( is.na(data.frame(datt)),1, sum)
NumberMissingByArray
# keep samples with missing values under 500
# in this case, all samples OK

# plots mean expression across all samples
pdf("sample mean expression.pdf",height=4, width=8)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:70), cex.names = 0.7)
dev.off()
# look for any obvious deviations in expression across samples

# sample dendrogram and trait heat map showing outliers
A=adjacency(t(datt),type="signed")                 #SELECT SIGNED OR UNSIGNED HERE
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(traits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(traits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
quartz()
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")
# the resulting plot shows a sample dendrogram and the spread of your traits across the cluster
# outlier samples will show as red in the outlierC row

# Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datt=datt[!remove.samples,]
traits=traits[!remove.samples,]
  
write.csv(traits, file="traits.csv")

save(datt,traits,file="wgcnaData.RData")


#### SOFT THRESHOLDS ####

library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
load("wgcnaData.RData")

# Try different betas ("soft threshold") - power factor for calling connections between genes
powers = c(seq(from = 2, to=26, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5,networkType="signed")

# Plot the results:
# Run from the line below to dev.off()
sizeGrWindow(9, 5)
pdf("soft threshold signed.pdf",height=4, width=8)

par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#### MAKING MODULES ####

# take a look at the threshold plots produced above, and the output table from the pickSoftThreshold command
# pick the power that corresponds with a SFT.R.sq value above 0.90
# if no powers result in value above 0.90, follow the table in the WGCNA FAQ: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

# run from the line below to the save command
s.th=12 # re-specify according to previous section
adjacency = adjacency(datt, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
# the grey module contains unassigned genes and is not considered a real module
# if you have error messages trying to generate the eigengene correlations, run this below
# check MEs, if grey shows NaN for all samples, then make sure to eliminate it using removeGreyME
# MEs = removeGreyME(MEs, greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules.RData")


#### MERGING MODULES ####

mm=load('1stPassModules.RData')
mm
lnames=load('wgcnaData.RData')
# traits
# head(datt)

quartz()

MEDissThres = 0 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
quartz()
# the grey module contains unassigned genes and is not considered a real module
# if you have error messages trying to generate the eigengene correlations, run this below
# check MEs, if grey shows NaN for all samples, then make sure to eliminate it using removeGreyME
# MEs = removeGreyME(MEs, greyMEName = paste(moduleColor.getMEprefix(), "grey", sep=""))
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)
# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed.RData")


#### MODULE CORRELATIONS ####
# plotting correlations with traits:
load(file = "networkdata_signed.RData")
load(file = "wgcnaData.RData");

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# module-trait correlations
quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

table(moduleColors) # gives numbers of genes in each module

# shows only significant correlations
quartz()
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor,2)
textMatrix = cors;
# paste(cors, "\n(",ps, ")", sep = "");
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-0.7,0.7),
               main = paste("Breviolum Module-Trait correlations"))

# module size barplot
labelShift=1000 # increase to move module size labels to the right
quartz()
par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,20000),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

# If it was first pass with no module merging, this is where you examine your heatmap and dendrogram of module eigengenes to see where you would like to set cut height (MEDissThres parameter) in the previous section to merge modules that are telling the same story for your trait data 
# A good way to do it is to find a group of similar modules in the heat map and then see at which tree height they connect in the dendrogram

#### NO MODULES ####
