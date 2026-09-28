#### packages ####

# install.packages("KOGMWU")
library(KOGMWU)

# loading KOG annotations
gene2kog=read.table("../../../Orbicella-faveolata-annotated-transcriptome/Ofaveolata_iso2kogClass.tab",sep="\t", fill=T)
head(gene2kog)


#### treatment all ####

LC_CC=load('../../outputs/transcriptomics/ofav/deseq2/LC_CC_lpv.RData')
LC_CC # names of datasets in the package
lpv.LC_CC=kog.mwu(LC_CC.p,gene2kog) 
lpv.LC_CC 

CH_CC=load('../../outputs/transcriptomics/ofav/deseq2/CH_CC_lpv.RData')
CH_CC # names of datasets in the package
lpv.CH_CC=kog.mwu(CH_CC.p,gene2kog) 
lpv.CH_CC

LH_CC=load('../../outputs/transcriptomics/ofav/deseq2/LH_CC_lpv.RData')
LH_CC # names of datasets in the package
lpv.LH_CC=kog.mwu(LH_CC.p,gene2kog) 
lpv.LH_CC 

CH_LC=load('../../outputs/transcriptomics/ofav/deseq2/CH_LC_lpv.RData')
CH_LC # names of datasets in the package
lpv.CH_LC=kog.mwu(CH_LC.p,gene2kog) 
lpv.CH_LC

LH_LC=load('../../outputs/transcriptomics/ofav/deseq2/LH_LC_lpv.RData')
LH_LC # names of datasets in the package
lpv.LH_LC=kog.mwu(LH_LC.p,gene2kog) 
lpv.LH_LC

LH_CH=load('../../outputs/transcriptomics/ofav/deseq2/LH_CH_lpv.RData')
LH_CH # names of datasets in the package
lpv.LH_CH=kog.mwu(LH_CH.p,gene2kog) 
lpv.LH_CH

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("LC_CC"=lpv.LC_CC,"CH_CC"=lpv.CH_CC,"LH_CC"=lpv.LH_CC,"CH_LC"=lpv.CH_LC,"LH_LC"=lpv.LH_LC,"LH_CH"=lpv.LH_CH))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white", filename="../../outputs/transcriptomics/ofav/kogmwu/KOG_urbanstress_ofav_host_treatment_lpv.pdf", width=7, height=8)

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="../../outputs/transcriptomics/ofav/kogmwu/KOG_urbanstress_ofav_host_treatment_corr_lpv.pdf", width=10, height=10)
par(mfrow=c(4,4))
corrPlot(x="LC_CC",y="CH_CC",ktable)
corrPlot(x="LC_CC",y="LH_CC",ktable)
corrPlot(x="LC_CC",y="CH_LC",ktable)
corrPlot(x="LC_CC",y="LH_LC",ktable)
corrPlot(x="LC_CC",y="LH_CH",ktable)
corrPlot(x="CH_CC",y="LH_CC",ktable)
corrPlot(x="CH_CC",y="CH_LC",ktable)
corrPlot(x="CH_CC",y="LH_LC",ktable)
corrPlot(x="CH_CC",y="LH_CH",ktable)
corrPlot(x="LH_CC",y="CH_LC",ktable)
corrPlot(x="LH_CC",y="LH_LC",ktable)
corrPlot(x="LH_CC",y="LH_CH",ktable)
corrPlot(x="CH_LC",y="LH_LC",ktable)
corrPlot(x="CH_LC",y="LH_CH",ktable)
corrPlot(x="LH_LC",y="LH_CH",ktable)
dev.off()


#### treatment filtered ####

LC_CC=load('../../outputs/transcriptomics/ofav/deseq2/LC_CC_lpv.RData')
LC_CC # names of datasets in the package
lpv.LC_CC=kog.mwu(LC_CC.p,gene2kog) 
lpv.LC_CC 

CH_CC=load('../../outputs/transcriptomics/ofav/deseq2/CH_CC_lpv.RData')
CH_CC # names of datasets in the package
lpv.CH_CC=kog.mwu(CH_CC.p,gene2kog) 
lpv.CH_CC

LH_CC=load('../../outputs/transcriptomics/ofav/deseq2/LH_CC_lpv.RData')
LH_CC # names of datasets in the package
lpv.LH_CC=kog.mwu(LH_CC.p,gene2kog) 
lpv.LH_CC 

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("LC_CC"=lpv.LC_CC,"CH_CC"=lpv.CH_CC,"LH_CC"=lpv.LH_CC))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white", filename="../../outputs/transcriptomics/ofav/kogmwu/KOG_urbanstress_ofav_host_treatment_filtered_lpv.pdf", width=7, height=8)

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="../../outputs/transcriptomics/ofav/kogmwu/KOG_urbanstress_ofav_host_treatment_filtered_corr_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="LC_CC",y="CH_CC",ktable)
corrPlot(x="LC_CC",y="LH_CC",ktable)
corrPlot(x="CH_CC",y="LH_CC",ktable)
dev.off()


#### urban vs reef ####

# urban vs reef DEGs come from the colony-level limma test in the DESeq2 script (urban_reef_lpv)
urban_reef=load('../../outputs/transcriptomics/ofav/deseq2/urban_reef_lpv.RData')
urban_reef # names of datasets in the package
lpv.urban_reef=kog.mwu(urban_reef.p,gene2kog)
lpv.urban_reef

# loading stress treatment results for comparison
LC_CC=load('../../outputs/transcriptomics/ofav/deseq2/LC_CC_lpv.RData')
lpv.LC_CC=kog.mwu(LC_CC.p,gene2kog)
CH_CC=load('../../outputs/transcriptomics/ofav/deseq2/CH_CC_lpv.RData')
lpv.CH_CC=kog.mwu(CH_CC.p,gene2kog)
LH_CC=load('../../outputs/transcriptomics/ofav/deseq2/LH_CC_lpv.RData')
lpv.LH_CC=kog.mwu(LH_CC.p,gene2kog)

# compiling a table of delta-ranks: does being urban resemble a stress response?
ktable=makeDeltaRanksTable(list("urban_reef"=lpv.urban_reef,"LC_CC"=lpv.LC_CC,"CH_CC"=lpv.CH_CC,"LH_CC"=lpv.LH_CC))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15, border_color="white", filename="../../outputs/transcriptomics/ofav/kogmwu/KOG_urbanstress_ofav_host_urban_reef_lpv.pdf", width=7, height=8)

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# creating a pub-ready corr plot
pdf(file="../../outputs/transcriptomics/ofav/kogmwu/KOG_urbanstress_ofav_host_urban_reef_corr_lpv.pdf", width=7, height=2.5)
par(mfrow=c(1,3))
corrPlot(x="urban_reef",y="LC_CC",ktable)
corrPlot(x="urban_reef",y="CH_CC",ktable)
corrPlot(x="urban_reef",y="LH_CC",ktable)
dev.off()
