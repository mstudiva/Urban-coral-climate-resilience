#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### DESEQ IMPORT TREATMENT ####

load("commongenes_DEGs.RData") # if previously run
ofav_LC_CC_lpv <- read.csv(file = "../../DESeq2/ofav/sym/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_CH_CC_lpv <- read.csv(file = "../../DESeq2/ofav/sym/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_CC_lpv <- read.csv(file = "../../DESeq2/ofav/sym/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_CH_LC_lpv <- read.csv(file = "../../DESeq2/ofav/sym/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_CH_lpv <- read.csv(file = "../../DESeq2/ofav/sym/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_LC_lpv <- read.csv(file = "../../DESeq2/ofav/sym/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ssid_LC_CC_lpv <- read.csv(file = "../../DESeq2/ssid/sym/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_CH_CC_lpv <- read.csv(file = "../../DESeq2/ssid/sym/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_CC_lpv <- read.csv(file = "../../DESeq2/ssid/sym/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_CH_LC_lpv <- read.csv(file = "../../DESeq2/ssid/sym/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_CH_lpv <- read.csv(file = "../../DESeq2/ssid/sym/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_LC_lpv <- read.csv(file = "../../DESeq2/ssid/sym/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)


#### DESEQ IMPORT SITE ####
ofav_Rainbow_Emerald_lpv <- read.csv(file = "../../DESeq2/ofav/sym/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_Star_Emerald_lpv <- read.csv(file = "../../DESeq2/ofav/sym/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Emerald_lpv <- read.csv(file = "../../DESeq2/ofav/sym/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_Star_Rainbow_lpv <- read.csv(file = "../../DESeq2/ofav/sym/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Rainbow_lpv <- read.csv(file = "../../DESeq2/ofav/sym/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Star_lpv <- read.csv(file = "../../DESeq2/ofav/sym/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ssid_Rainbow_Emerald_lpv <- read.csv(file = "../../DESeq2/ssid/sym/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_Star_Emerald_lpv <- read.csv(file = "../../DESeq2/ssid/sym/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Emerald_lpv <- read.csv(file = "../../DESeq2/ssid/sym/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_Star_Rainbow_lpv <- read.csv(file = "../../DESeq2/ssid/sym/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Rainbow_lpv <- read.csv(file = "../../DESeq2/ssid/sym/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Star_lpv <- read.csv(file = "../../DESeq2/ssid/sym/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)


#### DEG MATCHING TREATMENT ####

# This section of code does several things: 1) join -log10(pval) across species, 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3) adds gene annotations, and 4) then pulls on corresponding KOG classes

# LC vs CC for both species
ofav_LC_CC_lpv %>%
  inner_join(ssid_LC_CC_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="LC_CC", .before="gene") -> LC_CC

# CH vs CC for both species
ofav_CH_CC_lpv %>%
  inner_join(ssid_CH_CC_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="CH_CC", .before="gene") -> CH_CC

# LH vs CC for both species
ofav_LH_CC_lpv %>%
  inner_join(ssid_LH_CC_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="LH_CC", .before="gene") -> LH_CC

# CH vs LC for both species
ofav_CH_LC_lpv %>%
  inner_join(ssid_CH_LC_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="CH_LC", .before="gene") -> CH_LC

# LH vs CH for both species
ofav_LH_CH_lpv %>%
  inner_join(ssid_LH_CH_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="LH_CH", .before="gene") -> LH_CH

# LH vs LC for both species
ofav_LH_LC_lpv %>%
  inner_join(ssid_LH_LC_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="LH_LC", .before="gene") -> LH_LC

# joining all matching DEGs into a single dataframe
commongenes_treatment <- bind_rows(LC_CC,CH_CC,LH_CC,CH_LC,LH_CH,LH_LC)
write.csv(commongenes_treatment, file="commongenes_treatment.csv")


#### DEG MATCHING SITE ####

# This section of code does several things: 1) join -log10(pval) across species, 2) filter by 0.1 pval cutoff (log10(0.1)=1), 3) adds gene annotations, and 4) then pulls on corresponding KOG classes

# Rainbow vs Emerald for both species
ofav_Rainbow_Emerald_lpv %>%
  inner_join(ssid_Rainbow_Emerald_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="Rainbow_Emerald", .before="gene") -> Rainbow_Emerald

# Star vs Emerald for both species
ofav_Star_Emerald_lpv %>%
  inner_join(ssid_Star_Emerald_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="Star_Emerald", .before="gene") -> Star_Emerald

# MacN vs Emerald for both species
ofav_MacN_Emerald_lpv %>%
  inner_join(ssid_MacN_Emerald_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="MacN_Emerald", .before="gene") -> MacN_Emerald

# Star vs Rainbow for both species
ofav_Star_Rainbow_lpv %>%
  inner_join(ssid_Star_Rainbow_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="Star_Rainbow", .before="gene") -> Star_Rainbow

# MacN vs Star for both species
ofav_MacN_Star_lpv %>%
  inner_join(ssid_MacN_Star_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="MacN_Star", .before="gene") -> MacN_Star

# MacN vs Rainbow for both species
ofav_MacN_Rainbow_lpv %>%
  inner_join(ssid_MacN_Rainbow_lpv, by = c("gene" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>%
  mutate(comparison="MacN_Rainbow", .before="gene") -> MacN_Rainbow

# joining all matching DEGs into a single dataframe
commongenes_site <- bind_rows(Rainbow_Emerald,Star_Emerald,MacN_Emerald,Star_Rainbow,MacN_Rainbow,MacN_Star)
write.csv(commongenes_site, file="commongenes_site.csv")


#### VENN DIAGRAMS TREATMENT ####

# first creating a set of up/downregulated DEGs by species
CH_CC %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

CH_CC %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

CH_CC %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

CH_CC %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_CH_CC.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
LH_CC %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

LH_CC %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

LH_CC %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

LH_CC %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_LH_CC.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
CH_LC %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

CH_LC %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

CH_LC %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

CH_LC %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_CH_LC.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
LH_LC %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

LH_LC %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

LH_LC %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

LH_LC %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_LH_LC.pdf", height=10, width=12)
grid.draw(venn)
dev.off()


#### VENN DIAGRAMS SITE ####

# first creating a set of up/downregulated DEGs by species
Rainbow_Emerald %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

Rainbow_Emerald %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

Rainbow_Emerald %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

Rainbow_Emerald %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_Rainbow_Emerald.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
Star_Emerald %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

Star_Emerald %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

Star_Emerald %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

Star_Emerald %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_Star_Emerald.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
MacN_Emerald %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

MacN_Emerald %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

MacN_Emerald %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

MacN_Emerald %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_MacN_Emerald.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
Star_Rainbow %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

Star_Rainbow %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

Star_Rainbow %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

Star_Rainbow %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_Star_Rainbow.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
MacN_Rainbow %>%
  filter(lpv_ssid >= 1) %>%
  pull(gene) -> ssid_up

MacN_Rainbow %>%
  filter(lpv_ssid <= -1) %>%
  pull(gene) -> ssid_down

MacN_Rainbow %>%
  filter(lpv_ofav >= 1) %>%
  pull(gene) -> ofav_up

MacN_Rainbow %>%
  filter(lpv_ofav <= -1) %>%
  pull(gene) -> ofav_down

venn=venn.diagram(
  x = list("Ssid up"=ssid_up, "Ssid down"=ssid_down,"Ofav up"=ofav_up, "Ofav down"=ofav_down),
  filename=NULL,
  col = "transparent",
  fill = c("#ca0020", "#0571b0", "#f4a582", "#92c5de"),
  alpha = 0.5,
  label.col = c("red3","white","cornflowerblue","black","white","white","white", "black","darkred","grey25","white","white","grey25","darkblue","white"),
  cex = 3.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col =c("darkred", "darkblue", "red3", "cornflowerblue"),
  cat.cex = 3.5,
  cat.fontfamily = "sans",
  cat.just = list(c(0,0.5),c(0.75,0.5),c(0.5,0.5),c(0.5,0.5))
)
pdf(file="Venn_MacN_Rainbow.pdf", height=10, width=12)
grid.draw(venn)
dev.off()


#### HEATMAP DATA SUBSET ####

# first creating a column of combined gene names from both species, then removing unannotated genes
CH_CC %>%
  filter(!is.na(annot)) ->  CH_CC_heatmap

LH_CC %>%
  filter(!is.na(annot)) ->  LH_CC_heatmap

CH_LC %>%
  filter(!is.na(annot)) ->  CH_LC_heatmap

LH_LC %>%
  filter(!is.na(annot)) ->  LH_LC_heatmap

Rainbow_Emerald %>%
  filter(!is.na(annot)) ->  Rainbow_Emerald_heatmap

Star_Emerald %>%
  filter(!is.na(annot)) ->  Star_Emerald_heatmap

MacN_Emerald %>%
  filter(!is.na(annot)) ->  MacN_Emerald_heatmap

MacN_Rainbow %>%
  filter(!is.na(annot)) ->  MacN_Rainbow_heatmap


#### HEATMAPS TREATMENT ####

# only plotting comparisons with >15 shared, annotated DEGs
# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% CH_CC_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% CH_CC_heatmap$Protein_ofav, CH_CC_heatmap$gene, CH_CC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% CH_CC_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% CH_CC_heatmap$Protein_ssid, CH_CC_heatmap$gene, CH_CC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(CH_CC_heatmap$gene, CH_CC_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_CH_CC_p1e-6.pdf", height=10, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(CH_CC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% LH_CC_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% LH_CC_heatmap$Protein_ofav, LH_CC_heatmap$gene, LH_CC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% LH_CC_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% LH_CC_heatmap$Protein_ssid, LH_CC_heatmap$gene, LH_CC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(LH_CC_heatmap$gene, LH_CC_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_LH_CC_p1e-6.pdf", height=10, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(LH_CC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% CH_LC_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% CH_LC_heatmap$Protein_ofav, CH_LC_heatmap$gene, CH_LC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% CH_LC_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% CH_LC_heatmap$Protein_ssid, CH_LC_heatmap$gene, CH_LC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(CH_LC_heatmap$gene, CH_LC_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_CH_LC_p1e-6.pdf", height=10.5, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(CH_LC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% LH_LC_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% LH_LC_heatmap$Protein_ofav, LH_LC_heatmap$gene, LH_LC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% LH_LC_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% LH_LC_heatmap$Protein_ssid, LH_LC_heatmap$gene, LH_LC_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(LH_LC_heatmap$gene, LH_LC_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_LH_LC_p1e-6.pdf", height=10, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(LH_LC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()


#### HEATMAPS SITE ####

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% Rainbow_Emerald_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% Rainbow_Emerald_heatmap$Protein_ofav, Rainbow_Emerald_heatmap$gene, Rainbow_Emerald_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% Rainbow_Emerald_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% Rainbow_Emerald_heatmap$Protein_ssid, Rainbow_Emerald_heatmap$gene, Rainbow_Emerald_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(Rainbow_Emerald_heatmap$gene, Rainbow_Emerald_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_Rainbow_Emerald_p0.001.pdf", height=7, width=35)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(Rainbow_Emerald_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-3, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% Star_Emerald_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% Star_Emerald_heatmap$Protein_ofav, Star_Emerald_heatmap$gene, Star_Emerald_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% Star_Emerald_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% Star_Emerald_heatmap$Protein_ssid, Star_Emerald_heatmap$gene, Star_Emerald_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(Star_Emerald_heatmap$gene, Star_Emerald_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_Star_Emerald_p1e-6.pdf", height=5, width=37)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(Star_Emerald_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% MacN_Emerald_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% MacN_Emerald_heatmap$Protein_ofav, MacN_Emerald_heatmap$gene, MacN_Emerald_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% MacN_Emerald_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% MacN_Emerald_heatmap$Protein_ssid, MacN_Emerald_heatmap$gene, MacN_Emerald_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(MacN_Emerald_heatmap$gene, MacN_Emerald_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 1e-6
pdf(file="heatmap_MacN_Emerald_p1e-6.pdf", height=2.5, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(MacN_Emerald_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts
load("../../DESeq2/ofav/sym/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% MacN_Rainbow_heatmap$gene)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% MacN_Rainbow_heatmap$Protein_ofav, MacN_Rainbow_heatmap$gene, MacN_Rainbow_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/sym/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% MacN_Rainbow_heatmap$gene)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% MacN_Rainbow_heatmap$Protein_ssid, MacN_Rainbow_heatmap$gene, MacN_Rainbow_heatmap$gene)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of gene ID to gene annotations
gene_names <- as.data.frame(cbind(MacN_Rainbow_heatmap$gene, MacN_Rainbow_heatmap$annot))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_MacN_Rainbow_p0.1.pdf", height=2.25, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(MacN_Rainbow_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# saving dataframes
save(LC_CC, CH_CC, LH_CC, CH_LC, LH_CH, LH_LC, Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star, LC_CC_heatmap, CH_CC_heatmap, LH_CC_heatmap, CH_LC_heatmap, LH_LC_heatmap, Rainbow_Emerald_heatmap, Star_Emerald_heatmap, MacN_Emerald_heatmap, MacN_Rainbow_heatmap, file = "commongenes_DEGs.RData")


#### KOG MATCHING TREATMENT ####

load("commongenes_DEGs.RData") # if previously run

# no matching genes for LC vs CC in both species

# CH vs CC for both species
CH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

CH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

CH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

CH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_CH_CC.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# LH vs CC for both species
LH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

LH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

LH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

LH_CC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_LH_CC.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# CH vs LC for both species
CH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

CH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

CH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

CH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_CH_LC.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# LH vs LC for both species
LH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

LH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

LH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

LH_LC %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_LH_LC.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# no matching genes for LH vs CH in both species


#### KOG MATCHING SITE ####

# Rainbow vs Emerald for both species
Rainbow_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

Rainbow_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

Rainbow_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

Rainbow_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_Rainbow_Emerald.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# Star vs Emerald for both species
Star_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

Star_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

Star_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

Star_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_Star_Emerald.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# MacN vs Emerald for both species
MacN_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_up" = n) -> KOG_ofav_up

MacN_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ofav <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ofav_down" = n) -> KOG_ofav_down

MacN_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid >= 1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_up" = n) -> KOG_ssid_up

MacN_Emerald %>%
  mutate(KOG = replace(KOG, KOG == "", NA)) %>%
  filter(lpv_ssid <= -1) %>%
  count(KOG) %>%
  rename("KOG" = KOG, "ssid_down" = n) -> KOG_ssid_down

# joining all KOG class sums in a single dataframe
KOG_ofav_up %>%
  inner_join(KOG_ssid_up, by = "KOG") %>%
  inner_join(KOG_ofav_down, by = "KOG") %>%
  inner_join(KOG_ssid_down, by = "KOG") -> KOG_match

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
ggsave("KOG_MacN_Emerald.pdf", plot= KOG_sum, width=8, height=6, units="in", dpi=300)

# too few matching genes for Star vs Rainbow in both species

# no matching genes for MacN vs Rainbow in both species

# no matching genes for MacN vs Star in both species
