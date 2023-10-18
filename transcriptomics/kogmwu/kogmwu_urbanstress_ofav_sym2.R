#### packages ####

# install.packages("KOGMWU")
library(KOGMWU)

# loading KOG annotations
gene2kog=read.table("../../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",sep="\t", fill=T)
head(gene2kog)


#### treatment all ####

LC_CC=load('../../../DESeq2/ofav/sym2/LC_CC_lpv.RData')
LC_CC # names of datasets in the package
lpv.LC_CC=kog.mwu(LC_CC.p,gene2kog) 
lpv.LC_CC 

CH_CC=load('../../../DESeq2/ofav/sym2/CH_CC_lpv.RData')
CH_CC # names of datasets in the package
lpv.CH_CC=kog.mwu(CH_CC.p,gene2kog) 
lpv.CH_CC

LH_CC=load('../../../DESeq2/ofav/sym2/LH_CC_lpv.RData')
LH_CC # names of datasets in the package
lpv.LH_CC=kog.mwu(LH_CC.p,gene2kog) 
lpv.LH_CC 

CH_LC=load('../../../DESeq2/ofav/sym2/CH_LC_lpv.RData')
CH_LC # names of datasets in the package
lpv.CH_LC=kog.mwu(CH_LC.p,gene2kog) 
lpv.CH_LC

LH_LC=load('../../../DESeq2/ofav/sym2/LH_LC_lpv.RData')
LH_LC # names of datasets in the package
lpv.LH_LC=kog.mwu(LH_LC.p,gene2kog) 
lpv.LH_LC

LH_CH=load('../../../DESeq2/ofav/sym2/LH_CH_lpv.RData')
LH_CH # names of datasets in the package
lpv.LH_CH=kog.mwu(LH_CH.p,gene2kog) 
lpv.LH_CH

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("LC_CC"=lpv.LC_CC,"CH_CC"=lpv.CH_CC,"LH_CC"=lpv.LH_CC,"CH_LC"=lpv.CH_LC,"LH_LC"=lpv.LH_LC,"LH_CH"=lpv.LH_CH))
# many of the ranks are zero or NA for KOG classes, do not proceed with treatment comparisons


#### site all ####

Rainbow_Emerald=load('../../../DESeq2/ofav/sym2/Rainbow_Emerald_lpv.RData')
Rainbow_Emerald # names of datasets in the package
lpv.Rainbow_Emerald=kog.mwu(Rainbow_Emerald.p,gene2kog) 
lpv.Rainbow_Emerald 

Star_Emerald=load('../../../DESeq2/ofav/sym2/Star_Emerald_lpv.RData')
Star_Emerald # names of datasets in the package
lpv.Star_Emerald=kog.mwu(Star_Emerald.p,gene2kog) 
lpv.Star_Emerald 

MacN_Emerald=load('../../../DESeq2/ofav/sym2/MacN_Emerald_lpv.RData')
MacN_Emerald # names of datasets in the package
lpv.MacN_Emerald=kog.mwu(MacN_Emerald.p,gene2kog) 
lpv.MacN_Emerald 

Star_Rainbow=load('../../../DESeq2/ofav/sym2/Star_Rainbow_lpv.RData')
Star_Rainbow # names of datasets in the package
lpv.Star_Rainbow=kog.mwu(Star_Rainbow.p,gene2kog) 
lpv.Star_Rainbow 

MacN_Rainbow=load('../../../DESeq2/ofav/sym2/MacN_Rainbow_lpv.RData')
MacN_Rainbow # names of datasets in the package
lpv.MacN_Rainbow=kog.mwu(MacN_Rainbow.p,gene2kog) 
lpv.MacN_Rainbow 

MacN_Star=load('../../../DESeq2/ofav/sym2/MacN_Star_lpv.RData')
MacN_Star # names of datasets in the package
lpv.MacN_Star=kog.mwu(MacN_Star.p,gene2kog) 
lpv.MacN_Star 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Rainbow_Emerald"=lpv.Rainbow_Emerald,"Star_Emerald"=lpv.Star_Emerald,"MacN_Emerald"=lpv.MacN_Emerald,"Star_Rainbow"=lpv.Star_Rainbow,"MacN_Rainbow"=lpv.MacN_Rainbow,"MacN_Star"=lpv.MacN_Star))
ktable[is.na(ktable)] <- 0 # replacing any NAs in the rank table with zeroes

# Making a heatmap with hierarStarical clustering trees: 
pdf(file="KOG_urbanstress_ofav_sym2_site_lpv.pdf", width=7, height=7)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, celMacNeight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_urbanstress_ofav_sym2_site_corr_lpv.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="Rainbow_Emerald",y="Star_Emerald",ktable)
corrPlot(x="Rainbow_Emerald",y="MacN_Emerald",ktable)
corrPlot(x="Rainbow_Emerald",y="Star_Rainbow",ktable)
corrPlot(x="Rainbow_Emerald",y="MacN_Rainbow",ktable)
corrPlot(x="Rainbow_Emerald",y="MacN_Star",ktable)
corrPlot(x="Star_Emerald",y="MacN_Emerald",ktable)
corrPlot(x="Star_Emerald",y="Star_Rainbow",ktable)
corrPlot(x="Star_Emerald",y="MacN_Rainbow",ktable)
corrPlot(x="Star_Emerald",y="MacN_Star",ktable)
corrPlot(x="MacN_Emerald",y="Star_Rainbow",ktable)
corrPlot(x="MacN_Emerald",y="MacN_Rainbow",ktable)
corrPlot(x="MacN_Emerald",y="MacN_Star",ktable)
corrPlot(x="Star_Rainbow",y="MacN_Rainbow",ktable)
corrPlot(x="Star_Rainbow",y="MacN_Star",ktable)
corrPlot(x="MacN_Rainbow",y="MacN_Star",ktable)
dev.off()


#### treatment filtered ####

LC_CC=load('../../../DESeq2/ofav/sym2/LC_CC_lpv.RData')
LC_CC # names of datasets in the package
lpv.LC_CC=kog.mwu(LC_CC.p,gene2kog) 
lpv.LC_CC 

CH_CC=load('../../../DESeq2/ofav/sym2/CH_CC_lpv.RData')
CH_CC # names of datasets in the package
lpv.CH_CC=kog.mwu(CH_CC.p,gene2kog) 
lpv.CH_CC

LH_CC=load('../../../DESeq2/ofav/sym2/LH_CC_lpv.RData')
LH_CC # names of datasets in the package
lpv.LH_CC=kog.mwu(LH_CC.p,gene2kog) 
lpv.LH_CC 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("CH_CC"=lpv.CH_CC,"LH_CC"=lpv.LH_CC)) # removed LC_CC due to all ranks=0
ktable[is.na(ktable)] <- 0 # replacing any NAs in the rank table with zeroes

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pdf(file="KOG_urbanstress_ofav_sym2_treatment_filtered_lpv.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_urbanstress_ofav_sym2_treatment_filtered_corr_lpv.pdf", width=5, height=5)
par(mfrow=c(1,1))
corrPlot(x="CH_CC",y="LH_CC",ktable)
dev.off()


#### site filtered ####

Rainbow_Emerald=load('../../../DESeq2/ofav/sym2/Rainbow_Emerald_lpv.RData')
Rainbow_Emerald # names of datasets in the package
lpv.Rainbow_Emerald=kog.mwu(Rainbow_Emerald.p,gene2kog) 
lpv.Rainbow_Emerald 

Star_Emerald=load('../../../DESeq2/ofav/sym2/Star_Emerald_lpv.RData')
Star_Emerald # names of datasets in the package
lpv.Star_Emerald=kog.mwu(Star_Emerald.p,gene2kog) 
lpv.Star_Emerald 

MacN_Emerald=load('../../../DESeq2/ofav/sym2/MacN_Emerald_lpv.RData')
MacN_Emerald # names of datasets in the package
lpv.MacN_Emerald=kog.mwu(MacN_Emerald.p,gene2kog) 
lpv.MacN_Emerald 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Rainbow_Emerald"=lpv.Rainbow_Emerald,"Star_Emerald"=lpv.Star_Emerald,"MacN_Emerald"=lpv.MacN_Emerald))

# Making a heatmap with hierarStarical clustering trees: 
pdf(file="KOG_urbanstress_ofav_sym2_site_filtered_lpv.pdf", width=7, height=7)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, celMacNeight=15, border_color="white") 
while (!is.null(dev.list()))  dev.off()

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="KOG_urbanstress_ofav_sym2_site_filtered_corr_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="Rainbow_Emerald",y="Star_Emerald",ktable)
corrPlot(x="Rainbow_Emerald",y="MacN_Emerald",ktable)
corrPlot(x="Star_Emerald",y="MacN_Emerald",ktable)
dev.off()
