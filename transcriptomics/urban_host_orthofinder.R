#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(patchwork)   
library(cowplot)
library(stringr)


#### ORTHOFINDER ####

# Install orthofinder on your local machine using the tutorials (https://davidemms.github.io/menu/tutorials.html)

# Copy your translated protein fasta files (_out_PRO.fas) that you want to compare into a directory called 'orthofinder'
# If you have not already filtered by the longest contig per isogroup (by using fasta2SBH.pl during transcriptome annotation), follow step 7 of tutorial 2 above

# Run the following command in Terminal: 'orthofinder -f orthofinder/'
# Check the number of genes assigned to orthogroups (e.g., 'OrthoFinder assigned 39416 genes (89.6% of total) to 13781 orthogroups')
# Ideally, it should be >80%


#### ORTHOLOGS ####

# if loading in previously run data
load("orthofinder_DEGs.RData") 

# if starting from scratch
orthologs <- read.table(file = "orthofinder/OrthoFinder/Results_Apr10/Orthologues/Orthologues_Ofaveolata_out_PRO/Ofaveolata_out_PRO__v__Ssiderea_out_PRO.tsv", sep = "\t", header = TRUE, quote="", fill=FALSE)


#### DESEQ IMPORT TREATMENT ####

ofav_LC_CC_lpv <- read.csv(file = "../DESeq2/ofav/host/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_CH_CC_lpv <- read.csv(file = "../DESeq2/ofav/host/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_CC_lpv <- read.csv(file = "../DESeq2/ofav/host/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_CH_LC_lpv <- read.csv(file = "../DESeq2/ofav/host/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_CH_lpv <- read.csv(file = "../DESeq2/ofav/host/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_LC_lpv <- read.csv(file = "../DESeq2/ofav/host/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ssid_LC_CC_lpv <- read.csv(file = "../DESeq2/ssid/host/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_CH_CC_lpv <- read.csv(file = "../DESeq2/ssid/host/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_CC_lpv <- read.csv(file = "../DESeq2/ssid/host/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_CH_LC_lpv <- read.csv(file = "../DESeq2/ssid/host/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_CH_lpv <- read.csv(file = "../DESeq2/ssid/host/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_LC_lpv <- read.csv(file = "../DESeq2/ssid/host/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)


#### DESEQ IMPORT SITE ####
ofav_Rainbow_Emerald_lpv <- read.csv(file = "../DESeq2/ofav/host/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_Star_Emerald_lpv <- read.csv(file = "../DESeq2/ofav/host/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Emerald_lpv <- read.csv(file = "../DESeq2/ofav/host/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_Star_Rainbow_lpv <- read.csv(file = "../DESeq2/ofav/host/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Rainbow_lpv <- read.csv(file = "../DESeq2/ofav/host/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Star_lpv <- read.csv(file = "../DESeq2/ofav/host/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ssid_Rainbow_Emerald_lpv <- read.csv(file = "../DESeq2/ssid/host/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_Star_Emerald_lpv <- read.csv(file = "../DESeq2/ssid/host/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Emerald_lpv <- read.csv(file = "../DESeq2/ssid/host/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_Star_Rainbow_lpv <- read.csv(file = "../DESeq2/ssid/host/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Rainbow_lpv <- read.csv(file = "../DESeq2/ssid/host/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Star_lpv <- read.csv(file = "../DESeq2/ssid/host/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)


#### DEG MATCHING TREATMENT ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds ofav and ssid gene annotations, and 5) then pulls on corresponding KOG classes

# LC vs CC for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LC_CC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LC_CC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="LC_CC", .before="Orthogroup") -> LC_CC
LC_CC$Orthogroup <- make.unique(LC_CC$Orthogroup, sep = "_") 

# CH vs CC for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_CH_CC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_CH_CC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="CH_CC", .before="Orthogroup") -> CH_CC
CH_CC$Orthogroup <- make.unique(CH_CC$Orthogroup, sep = "_") 

# LH vs CC for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LH_CC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LH_CC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="LH_CC", .before="Orthogroup") -> LH_CC
LH_CC$Orthogroup <- make.unique(LH_CC$Orthogroup, sep = "_") 

# CH vs LC for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_CH_LC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_CH_LC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="CH_LC", .before="Orthogroup") -> CH_LC
CH_LC$Orthogroup <- make.unique(CH_LC$Orthogroup, sep = "_") 

# LH vs CH for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LH_CH_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LH_CH_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="LH_CH", .before="Orthogroup") -> LH_CH
LH_CH$Orthogroup <- make.unique(LH_CH$Orthogroup, sep = "_") 

# LH vs LC for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LH_LC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LH_LC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="LH_LC", .before="Orthogroup") -> LH_LC
LH_LC$Orthogroup <- make.unique(LH_LC$Orthogroup, sep = "_") 

# joining all matching DEGs into a single dataframe
orthofinder_treatment <- bind_rows(LC_CC,CH_CC,LH_CC,CH_LC,LH_CH,LH_LC)
write.csv(orthofinder_treatment, file="orthofinder_treatment.csv")


#### DEG MATCHING SITE ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds ofav and ssid gene annotations, and 5) then pulls on corresponding KOG classes

# Rainbow vs Emerald for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_Rainbow_Emerald_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_Rainbow_Emerald_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="Rainbow_Emerald", .before="Orthogroup") -> Rainbow_Emerald
Rainbow_Emerald$Orthogroup <- make.unique(Rainbow_Emerald$Orthogroup, sep = "_") 

# Star vs Emerald for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_Star_Emerald_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_Star_Emerald_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="Star_Emerald", .before="Orthogroup") -> Star_Emerald
Star_Emerald$Orthogroup <- make.unique(Star_Emerald$Orthogroup, sep = "_") 

# MacN vs Emerald for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_MacN_Emerald_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_MacN_Emerald_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="MacN_Emerald", .before="Orthogroup") -> MacN_Emerald
MacN_Emerald$Orthogroup <- make.unique(MacN_Emerald$Orthogroup, sep = "_") 

# Star vs Rainbow for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_Star_Rainbow_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_Star_Rainbow_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="Star_Rainbow", .before="Orthogroup") -> Star_Rainbow
Star_Rainbow$Orthogroup <- make.unique(Star_Rainbow$Orthogroup, sep = "_") 

# MacN vs Rainbow for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_MacN_Rainbow_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_MacN_Rainbow_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="MacN_Rainbow", .before="Orthogroup") -> MacN_Rainbow
MacN_Rainbow$Orthogroup <- make.unique(MacN_Rainbow$Orthogroup, sep = "_") 

# MacN vs Star for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Ssiderea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_MacN_Star_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_MacN_Star_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../Annotations/ssid/locatelli/Ssiderea_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  mutate(comparison="MacN_Star", .before="Orthogroup") -> MacN_Star
MacN_Star$Orthogroup <- make.unique(MacN_Star$Orthogroup, sep = "_") 

# joining all matching DEGs into a single dataframe
orthofinder_site <- bind_rows(Rainbow_Emerald,Star_Emerald,MacN_Emerald,Star_Rainbow,MacN_Rainbow,MacN_Star)
write.csv(orthofinder_site, file="orthofinder_site.csv")


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

# --- Helpers ---
sgn <- function(z) ifelse(z > 0, 1, ifelse(z < 0, -1, 0))

clean_lpv_df <- function(df, cap = inf_cap) {
  d <- df
  d$lpv_ofav <- as.numeric(d$lpv_ofav)
  d$lpv_ssid <- as.numeric(d$lpv_ssid)
  # Cap infinities from p=0
  d$lpv_ofav[ is.infinite(d$lpv_ofav) & d$lpv_ofav > 0 ] <- cap
  d$lpv_ofav[ is.infinite(d$lpv_ofav) & d$lpv_ofav < 0 ] <- -cap
  d$lpv_ssid[ is.infinite(d$lpv_ssid) & d$lpv_ssid > 0 ] <- cap
  d$lpv_ssid[ is.infinite(d$lpv_ssid) & d$lpv_ssid < 0 ] <- -cap
  # Drop non-finite rows
  d <- d[is.finite(d$lpv_ofav) & is.finite(d$lpv_ssid), , drop = FALSE]
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

sig_category4 <- function(lpv_ofav, lpv_ssid, c1, c2, c3, levels_vec = SIG_LEVELS) {
  v <- sort(c(c1, c2, c3)); c1 <- v[1]; c2 <- v[2]; c3 <- v[3]
  lv <- pmin(abs(lpv_ofav), abs(lpv_ssid))  # min significance across species in |−log10 p|
  out <- cut(
    lv,
    breaks = c(-Inf, c1, c2, c3, Inf),
    labels = pvalue_labels(c1, c2, c3),
    include.lowest = TRUE, right = TRUE
  )
  factor(out, levels = levels_vec)
}

compute_global_limits <- function(comparisons, min_span = 0.6, pad_frac = 0.05) {
  vals <- unlist(lapply(comparisons, function(df) {
    d <- clean_lpv_df(df)
    c(d$lpv_ofav, d$lpv_ssid)
  }))
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs)) max_abs <- min_span / 2
  lim <- max(max_abs * (1 + pad_frac), min_span / 2)
  c(-lim, lim)
}

# --- One-panel constructor (returns a ggplot object) ---
make_panel <- function(df, title,
                       c1 = sig_cut1_global, c2 = sig_cut2_global, c3 = sig_cut3_global) {
  dat <- clean_lpv_df(df)   # Apply lpv threshold of 1 to plotted points (max across species)   dat <- subset(dat, pmax(abs(lpv_ofav), abs(lpv_ssid)) >= 1.0)
  # Apply lpv threshold of 1 to plotted points (max across species)
  dat <- subset(dat, pmax(abs(lpv_ofav), abs(lpv_ssid)) >= 1.0)
  
  # Relationship class
  same <- (dat$lpv_ofav >= 0 & dat$lpv_ssid >= 0) | (dat$lpv_ofav <= 0 & dat$lpv_ssid <= 0)
  dat$rel_class <- factor(ifelse(same, "Direct", "Inverse"), levels = REL_LEVELS)
  
  # Significance tiers
  dat$sig_cat <- sig_category4(dat$lpv_ofav, dat$lpv_ssid, c1, c2, c3, levels_vec = SIG_LEVELS)
  
  # --- Data-driven square limits (not forced around 0), with min span ---
  x_min <- min(dat$lpv_ofav, na.rm = TRUE); x_max <- max(dat$lpv_ofav, na.rm = TRUE)
  y_min <- min(dat$lpv_ssid, na.rm = TRUE); y_max <- max(dat$lpv_ssid, na.rm = TRUE)
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
  st <- safe_stats(dat$lpv_ofav, dat$lpv_ssid)
  fmt_num <- function(x) if (is.na(x)) "NA" else if (abs(x) < 1e-3) format(x, digits = 2, scientific = TRUE) else sprintf("%.3f", x)
  lab_text <- if (is.na(st$r)) sprintf("n = %d", st$n) else sprintf("r = %.2f, p = %s, n = %d", st$r, fmt_num(st$p), st$n)
  
  # Legend labels for line types
  ref_levels <- c("Direct 1:1", "Inverse 1:1")
  
  ggplot(dat, aes(x = lpv_ofav, y = lpv_ssid)) +
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
      x = paste0("O. faveolata (signed ", MINUS, "log10 p)"),
      y = paste0("S. siderea (signed ", MINUS, "log10 p)"),
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

# --- Build all six panels ---
# Expect these data frames to exist in your environment:
# LC_CC, CH_CC, LH_CC, CH_LC, LH_LC, LH_CH
comparisons <- list(
  "OA" = LC_CC,
  "Bleaching" = CH_CC,
  "OA + Bleaching" = LH_CC,
  "Bleaching vs OA" = CH_LC,
  "OA + Bleaching vs OA" = LH_LC,
  "OA + Bleaching vs Bleaching" = LH_CH
)

panels <- mapply(
  FUN = function(title, df) make_panel(df, title, sig_cut1_global, sig_cut2_global, sig_cut3_global),
  title = names(comparisons),
  df    = comparisons,
  SIMPLIFY = FALSE
)

# --- Hide redundant axis TITLES only (keep ticks & numbers) for 2x3 layout ---
for (i in seq_along(panels)) {
  row <- ceiling(i / 3)
  col <- i - (row - 1) * 3
  if (col != 1) panels[[i]] <- panels[[i]] + theme(axis.title.y = element_blank())
  if (row != 2) panels[[i]] <- panels[[i]] + theme(axis.title.x = element_blank())
}

# --- Single legend workflow ---
one_legend  <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
legend_plot <- cowplot::ggdraw(one_legend)

panels_noleg <- lapply(panels, function(p) p + theme(legend.position = "none"))

grid  <- wrap_plots(panels_noleg, ncol = 3)
panel <- grid | legend_plot
panel <- panel + plot_layout(widths = c(1, 0.30))

# --- Show & Save ---
print(panel)
# Use Cairo devices if you want true Unicode glyphs in output:
# ggsave("orthogroup_correlation_treatment.pdf", panel, width = 14, height = 9, device = cairo_pdf)
ggsave("orthogroup correlation treatment.pdf", panel, width = 14, height = 9)


#### CORRELATION DATAFRAME TREATMENT ####

# knobs
delta_cut_dd <- 2.0   # |Δ −log10 p| >= 2  => ≥100x p-value difference
ratio_cut_dd <- 100   # p-ratio >= 100
high_cut_dd  <- sig_cut2_global   # "strong" tier (default: 2.0 ~ p <= 0.01)
low_cut_dd   <- sig_cut1_global   # "weak"   tier (default: 1.3 ~ p <= 0.05)
min_sig      <- 1               # require at least this |−log10 p| in either species (0 = none)

classify_for_filters <- function(df) {
  d <- clean_lpv_df(df)
  same <- (d$lpv_ofav >= 0 & d$lpv_ssid >= 0) | (d$lpv_ofav <= 0 & d$lpv_ssid <= 0)
  d$rel_class      <- ifelse(same, "Direct", "Inverse")
  d$abs_ofav       <- abs(d$lpv_ofav)
  d$abs_ssid       <- abs(d$lpv_ssid)
  d$delta_log10p   <- abs(d$abs_ofav - d$abs_ssid)            # |Δ −log10 p|
  d$p_ofav         <- 10^(-d$abs_ofav)
  d$p_ssid         <- 10^(-d$abs_ssid)
  d$p_ratio        <- pmax(d$p_ofav, d$p_ssid) / pmin(d$p_ofav, d$p_ssid)  # >= 1
  d
}

# helper: align columns across data frames then rbind (base R only)
rbind_align <- function(...) {
  dfs <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(x) {
    if (is.null(x) || !nrow(x)) {
      # create empty with all cols
      y <- as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
      return(y)
    }
    miss <- setdiff(cols, names(x))
    for (m in miss) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, dfs2)
}

filter_one <- function(df, treatment_name) {
  d <- classify_for_filters(df)
  meets_min <- (pmax(d$abs_ofav, d$abs_ssid) >= min_sig)
  
  # 1) ALL inverse
  inverse <- d[d$rel_class == "Inverse" & meets_min, , drop = FALSE]
  if (nrow(inverse)) {
    inverse$filter_type           <- "Inverse"
    inverse$relationship_category <- "Inverse"
    inverse$discord_direction     <- NA_character_
    inverse$delta_cut_used        <- NA_real_
    inverse$high_cut_used         <- NA_real_
    inverse$low_cut_used          <- NA_real_
    inverse$ratio_cut_used        <- NA_real_
    inverse$treatment             <- treatment_name
  }
  
  # 2) STRICT direct-discordant
  strong_ofav_weak_ssid <- d$abs_ofav >= high_cut_dd & d$abs_ssid <= low_cut_dd
  strong_ssid_weak_ofav <- d$abs_ssid >= high_cut_dd & d$abs_ofav <= low_cut_dd
  cross_tier <- strong_ofav_weak_ssid | strong_ssid_weak_ofav
  
  dd_idx <- d$rel_class == "Direct" &
    d$delta_log10p >= delta_cut_dd &
    d$p_ratio      >= ratio_cut_dd &
    cross_tier &
    meets_min
  
  direct_discordant <- d[dd_idx, , drop = FALSE]
  if (nrow(direct_discordant)) {
    direct_discordant$filter_type           <- "DirectDiscordant"
    direct_discordant$relationship_category <- "Direct (discordant significance)"
    # compute direction row-wise inside the subset
    so <- direct_discordant$abs_ofav >= high_cut_dd & direct_discordant$abs_ssid <= low_cut_dd
    direct_discordant$discord_direction <- ifelse(so, "Ofav >> Ssid", "Ssid >> Ofav")
    direct_discordant$delta_cut_used    <- delta_cut_dd
    direct_discordant$high_cut_used     <- high_cut_dd
    direct_discordant$low_cut_used      <- low_cut_dd
    direct_discordant$ratio_cut_used    <- ratio_cut_dd
    direct_discordant$treatment         <- treatment_name
  }
  
  # align columns before binding
  rbind_align(inverse, direct_discordant)
}

# build ONE master dataframe
treatment_filtered_genes <- do.call(
  rbind,
  mapply(
    FUN = function(nm, df) filter_one(df, nm),
    nm  = names(comparisons),
    df  = comparisons,
    SIMPLIFY = FALSE
  )
)

# sanity checks
if (exists("treatment_filtered_genes") && nrow(treatment_filtered_genes) > 0) {
  cat("Rows total:", nrow(treatment_filtered_genes), "\n")
  print(table(treatment_filtered_genes$relationship_category, useNA = "ifany"))
  print(table(treatment_filtered_genes$treatment,   useNA = "ifany"))
} else {
  warning("treatment_filtered_genes is empty with current strict thresholds.")
}

# total counts
dplyr::count(treatment_filtered_genes, treatment, relationship_category)

write.csv(treatment_filtered_genes, file = "orthogroup correlation discordant treatment.csv")


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

# --- Helpers ---
sgn <- function(z) ifelse(z > 0, 1, ifelse(z < 0, -1, 0))

clean_lpv_df <- function(df, cap = inf_cap) {
  d <- df
  d$lpv_ofav <- as.numeric(d$lpv_ofav)
  d$lpv_ssid <- as.numeric(d$lpv_ssid)
  # Cap infinities from p=0
  d$lpv_ofav[ is.infinite(d$lpv_ofav) & d$lpv_ofav > 0 ] <- cap
  d$lpv_ofav[ is.infinite(d$lpv_ofav) & d$lpv_ofav < 0 ] <- -cap
  d$lpv_ssid[ is.infinite(d$lpv_ssid) & d$lpv_ssid > 0 ] <- cap
  d$lpv_ssid[ is.infinite(d$lpv_ssid) & d$lpv_ssid < 0 ] <- -cap
  # Drop non-finite rows
  d <- d[is.finite(d$lpv_ofav) & is.finite(d$lpv_ssid), , drop = FALSE]
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

sig_category4 <- function(lpv_ofav, lpv_ssid, c1, c2, c3, levels_vec = SIG_LEVELS) {
  v <- sort(c(c1, c2, c3)); c1 <- v[1]; c2 <- v[2]; c3 <- v[3]
  lv <- pmin(abs(lpv_ofav), abs(lpv_ssid))  # min significance across species in |−log10 p|
  out <- cut(
    lv,
    breaks = c(-Inf, c1, c2, c3, Inf),
    labels = pvalue_labels(c1, c2, c3),
    include.lowest = TRUE, right = TRUE
  )
  factor(out, levels = levels_vec)
}

compute_global_limits <- function(comparisons, min_span = 0.6, pad_frac = 0.05) {
  vals <- unlist(lapply(comparisons, function(df) {
    d <- clean_lpv_df(df)
    c(d$lpv_ofav, d$lpv_ssid)
  }))
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs)) max_abs <- min_span / 2
  lim <- max(max_abs * (1 + pad_frac), min_span / 2)
  c(-lim, lim)
}

# --- One-panel constructor (returns a ggplot object) ---
make_panel <- function(df, title,
                       c1 = sig_cut1_global, c2 = sig_cut2_global, c3 = sig_cut3_global) {
  dat <- clean_lpv_df(df)   # Apply lpv threshold of 1 to plotted points (max across species)   dat <- subset(dat, pmax(abs(lpv_ofav), abs(lpv_ssid)) >= 1.0)
  
  # Relationship class
  same <- (dat$lpv_ofav >= 0 & dat$lpv_ssid >= 0) | (dat$lpv_ofav <= 0 & dat$lpv_ssid <= 0)
  dat$rel_class <- factor(ifelse(same, "Direct", "Inverse"), levels = REL_LEVELS)
  
  # Significance tiers
  dat$sig_cat <- sig_category4(dat$lpv_ofav, dat$lpv_ssid, c1, c2, c3, levels_vec = SIG_LEVELS)
  
  # --- Data-driven square limits (not forced around 0), with min span ---
  x_min <- min(dat$lpv_ofav, na.rm = TRUE); x_max <- max(dat$lpv_ofav, na.rm = TRUE)
  y_min <- min(dat$lpv_ssid, na.rm = TRUE); y_max <- max(dat$lpv_ssid, na.rm = TRUE)
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
  st <- safe_stats(dat$lpv_ofav, dat$lpv_ssid)
  fmt_num <- function(x) if (is.na(x)) "NA" else if (abs(x) < 1e-3) format(x, digits = 2, scientific = TRUE) else sprintf("%.3f", x)
  lab_text <- if (is.na(st$r)) sprintf("n = %d", st$n) else sprintf("r = %.2f, p = %s, n = %d", st$r, fmt_num(st$p), st$n)
  
  # Legend labels for line types
  ref_levels <- c("Direct 1:1", "Inverse 1:1")
  
  ggplot(dat, aes(x = lpv_ofav, y = lpv_ssid)) +
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
      x = paste0("O. faveolata (signed ", MINUS, "log10 p)"),
      y = paste0("S. siderea (signed ", MINUS, "log10 p)"),
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

# --- Build all six panels ---
# Expect these data frames to exist in your environment:
# Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star
comparisons <- list(
  "MacArthur North vs Emerald Reef" = MacN_Emerald,
  "MacArthur North vs Rainbow Reef" = MacN_Rainbow,
  "MacArthur North vs Star Island" = MacN_Star,
  "Star Island vs Emerald Reef" = Star_Emerald,
  "Star Island vs Rainbow Reef" = Star_Rainbow,
  "Rainbow Reef vs Emerald Reef" = Rainbow_Emerald
)

panels <- mapply(
  FUN = function(title, df) make_panel(df, title, sig_cut1_global, sig_cut2_global, sig_cut3_global),
  title = names(comparisons),
  df    = comparisons,
  SIMPLIFY = FALSE
)

# --- Hide redundant axis TITLES only (keep ticks & numbers) for 2x3 layout ---
for (i in seq_along(panels)) {
  row <- ceiling(i / 3)
  col <- i - (row - 1) * 3
  if (col != 1) panels[[i]] <- panels[[i]] + theme(axis.title.y = element_blank())
  if (row != 2) panels[[i]] <- panels[[i]] + theme(axis.title.x = element_blank())
}

# --- Single legend workflow ---
one_legend  <- cowplot::get_legend(panels[[1]] + theme(legend.position = "right"))
legend_plot <- cowplot::ggdraw(one_legend)

panels_noleg <- lapply(panels, function(p) p + theme(legend.position = "none"))

grid  <- wrap_plots(panels_noleg, ncol = 3)
panel <- grid | legend_plot
panel <- panel + plot_layout(widths = c(1, 0.30))

# --- Show & Save ---
print(panel)
# Use Cairo devices if you want true Unicode glyphs in output:
# ggsave("orthogroup_correlation_site.pdf", panel, width = 14, height = 9, device = cairo_pdf)
ggsave("orthogroup correlation site.pdf", panel, width = 14, height = 9)


#### CORRELATION DATAFRAME SITE ####

# knobs
delta_cut_dd <- 2.0   # |Δ −log10 p| >= 2  => ≥100x p-value difference
ratio_cut_dd <- 100   # p-ratio >= 100
high_cut_dd  <- sig_cut2_global   # "strong" tier (default: 2.0 ~ p <= 0.01)
low_cut_dd   <- sig_cut1_global   # "weak"   tier (default: 1.3 ~ p <= 0.05)
min_sig      <- 1               # require at least this |−log10 p| in either species (0 = none)

classify_for_filters <- function(df) {
  d <- clean_lpv_df(df)
  same <- (d$lpv_ofav >= 0 & d$lpv_ssid >= 0) | (d$lpv_ofav <= 0 & d$lpv_ssid <= 0)
  d$rel_class      <- ifelse(same, "Direct", "Inverse")
  d$abs_ofav       <- abs(d$lpv_ofav)
  d$abs_ssid       <- abs(d$lpv_ssid)
  d$delta_log10p   <- abs(d$abs_ofav - d$abs_ssid)            # |Δ −log10 p|
  d$p_ofav         <- 10^(-d$abs_ofav)
  d$p_ssid         <- 10^(-d$abs_ssid)
  d$p_ratio        <- pmax(d$p_ofav, d$p_ssid) / pmin(d$p_ofav, d$p_ssid)  # >= 1
  d
}

# helper: align columns across data frames then rbind (base R only)
rbind_align <- function(...) {
  dfs <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(x) {
    if (is.null(x) || !nrow(x)) {
      # create empty with all cols
      y <- as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
      return(y)
    }
    miss <- setdiff(cols, names(x))
    for (m in miss) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, dfs2)
}

filter_one <- function(df, site_name) {
  d <- classify_for_filters(df)
  meets_min <- (pmax(d$abs_ofav, d$abs_ssid) >= min_sig)
  
  # 1) ALL inverse
  inverse <- d[d$rel_class == "Inverse" & meets_min, , drop = FALSE]
  if (nrow(inverse)) {
    inverse$filter_type           <- "Inverse"
    inverse$relationship_category <- "Inverse"
    inverse$discord_direction     <- NA_character_
    inverse$delta_cut_used        <- NA_real_
    inverse$high_cut_used         <- NA_real_
    inverse$low_cut_used          <- NA_real_
    inverse$ratio_cut_used        <- NA_real_
    inverse$site             <- site_name
  }
  
  # 2) STRICT direct-discordant
  strong_ofav_weak_ssid <- d$abs_ofav >= high_cut_dd & d$abs_ssid <= low_cut_dd
  strong_ssid_weak_ofav <- d$abs_ssid >= high_cut_dd & d$abs_ofav <= low_cut_dd
  cross_tier <- strong_ofav_weak_ssid | strong_ssid_weak_ofav
  
  dd_idx <- d$rel_class == "Direct" &
    d$delta_log10p >= delta_cut_dd &
    d$p_ratio      >= ratio_cut_dd &
    cross_tier &
    meets_min
  
  direct_discordant <- d[dd_idx, , drop = FALSE]
  if (nrow(direct_discordant)) {
    direct_discordant$filter_type           <- "DirectDiscordant"
    direct_discordant$relationship_category <- "Direct (discordant significance)"
    # compute direction row-wise inside the subset
    so <- direct_discordant$abs_ofav >= high_cut_dd & direct_discordant$abs_ssid <= low_cut_dd
    direct_discordant$discord_direction <- ifelse(so, "Ofav >> Ssid", "Ssid >> Ofav")
    direct_discordant$delta_cut_used    <- delta_cut_dd
    direct_discordant$high_cut_used     <- high_cut_dd
    direct_discordant$low_cut_used      <- low_cut_dd
    direct_discordant$ratio_cut_used    <- ratio_cut_dd
    direct_discordant$site         <- site_name
  }
  
  # align columns before binding
  rbind_align(inverse, direct_discordant)
}

# build ONE master dataframe
site_filtered_genes <- do.call(
  rbind,
  mapply(
    FUN = function(nm, df) filter_one(df, nm),
    nm  = names(comparisons),
    df  = comparisons,
    SIMPLIFY = FALSE
  )
)

# sanity checks
if (exists("site_filtered_genes") && nrow(site_filtered_genes) > 0) {
  cat("Rows total:", nrow(site_filtered_genes), "\n")
  print(table(site_filtered_genes$relationship_category, useNA = "ifany"))
  print(table(site_filtered_genes$site,   useNA = "ifany"))
} else {
  warning("site_filtered_genes is empty with current strict thresholds.")
}

# total counts
dplyr::count(site_filtered_genes, site, relationship_category)

write.csv(site_filtered_genes, file = "orthogroup correlation discordant site.csv")


#### SAVING DATAFRAMES ####

# saving dataframes
save(orthologs, LC_CC, CH_CC, LH_CC, CH_LC, LH_CH, LH_LC, Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star, file = "orthofinder_DEGs.RData")
