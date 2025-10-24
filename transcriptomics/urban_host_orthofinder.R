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


#### SAVING DATAFRAMES ####

# saving dataframes
save(orthologs, LC_CC, CH_CC, LH_CC, CH_LC, LH_CH, LH_LC, Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star, file = "orthofinder_DEGs.RData")


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


#### GENE BOXPLOTS ####

# ---- (tiny preamble) ----
while (dev.cur() > 1) try(dev.off(), silent = TRUE)
options(device = "pdf")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(DESeq2)
  library(SummarizedExperiment)
  library(readr)
  library(cowplot)
  library(grid)
  library(scales)
  library(rstatix)
  library(stringr)
  library(rlang)
})
`%||%` <- rlang::`%||%`

# --- Ensure dplyr functions override Bioconductor ones ---
for (fn in c("filter","select","mutate","arrange","summarize","slice"))
  assign(fn, get(fn, envir = asNamespace("dplyr")), envir = globalenv())

set.seed(42)
options(ggplot2.useDingbats = FALSE)

# ---- Helpers ----
.safe <- function(x) gsub("[^[:alnum:]_\\-\\.]+", "_", x)

# ---- Top-height calculator for brackets ----
.bracket_required_top <- function(df, tukey_tbl, lvls,
                                  pad_dex   = 0.15,
                                  gap_dex   = 0.15,
                                  lbl_dex   = 0.005,
                                  min_sep_dex = 0.08,
                                  micro_dex = 0.020) {
  if (is.null(tukey_tbl) || !nrow(tukey_tbl)) return(NA_real_)
  tuk <- tukey_tbl
  tuk$group1 <- as.character(tuk$group1)
  tuk$group2 <- as.character(tuk$group2)
  tuk <- tuk[tuk$group1 %in% lvls & tuk$group2 %in% lvls, , drop=FALSE]
  if (!nrow(tuk)) return(NA_real_)
  if ("p.adj.signif" %in% names(tuk)) tuk <- dplyr::filter(tuk, !is.na(p.adj.signif) & tolower(p.adj.signif) != "ns")
  else if ("p.adj" %in% names(tuk))   tuk <- dplyr::filter(tuk, !is.na(p.adj) & p.adj < 0.05)
  else if ("p" %in% names(tuk))       tuk <- dplyr::filter(tuk, !is.na(p) & p < 0.05)
  else return(NA_real_)
  if (!nrow(tuk)) return(NA_real_)
  
  x1 <- match(tuk$group1, lvls); x2 <- match(tuk$group2, lvls)
  span <- abs(x2 - x1)
  tuk  <- tuk[order(span), , drop=FALSE]
  
  overlaps <- function(a1,a2,b1,b2) !(a2 < b1 || b2 < a1)
  lanes <- list(); lane_id <- integer(nrow(tuk))
  for (i in seq_len(nrow(tuk))) {
    a1 <- min(x1[i], x2[i]); a2 <- max(x1[i], x2[i]); placed <- FALSE
    for (L in seq_along(lanes)) {
      lane_ok <- TRUE
      for (j in seq_len(nrow(lanes[[L]]))) if (overlaps(a1,a2, lanes[[L]]$a1[j], lanes[[L]]$a2[j])) { lane_ok <- FALSE; break }
      if (lane_ok) { lanes[[L]] <- rbind(lanes[[L]], data.frame(a1=a1,a2=a2)); lane_id[i] <- L; placed <- TRUE; break }
    }
    if (!placed) { lanes[[length(lanes)+1]] <- data.frame(a1=a1,a2=a2); lane_id[i] <- length(lanes) }
  }
  
  df_pos <- df %>% dplyr::filter(is.finite(count), count > 0)
  group_max <- tapply(df_pos$count, factor(df_pos$grp, levels=lvls), max, na.rm=TRUE)
  group_max[is.infinite(group_max) | is.na(group_max)] <- 1
  base_h   <- pmax(group_max[tuk$group1], group_max[tuk$group2])
  base_log <- log10(pmax(base_h, 1))
  
  within_lane_idx <- ave(seq_along(lane_id), lane_id, FUN=function(ix) rank(ix, ties.method="first"))
  y_log   <- base_log + pad_dex + (lane_id - 1)*gap_dex + (within_lane_idx - 1)*pmax(micro_dex, min_sep_dex)
  lbl_log <- y_log + lbl_dex
  10^(max(lbl_log, na.rm=TRUE) + 0.010)
}

# ---- Draw significant Tukey brackets ----
.add_sig_brackets <- function(p, tukey_tbl, lvls, top_y, df,
                              pad_dex=0.2, gap_dex=0.15, tip_dex=0.025,
                              lbl_dex=0.005, cap_dex=0.010, min_clr=0.06, min_sep_dex=0.08){
  if (is.null(tukey_tbl) || !nrow(tukey_tbl)) return(p)
  tuk <- tukey_tbl
  tuk$group1 <- as.character(tuk$group1); tuk$group2 <- as.character(tuk$group2)
  tuk <- tuk[tuk$group1 %in% lvls & tuk$group2 %in% lvls, , drop = FALSE]
  if (!nrow(tuk)) return(p)
  if ("p.adj.signif" %in% names(tuk)) tuk <- dplyr::filter(tuk, !is.na(p.adj.signif) & tolower(p.adj.signif) != "ns")
  else if ("p.adj" %in% names(tuk))   tuk <- dplyr::filter(tuk, !is.na(p.adj) & p.adj < 0.05)
  else if ("p" %in% names(tuk))       tuk <- dplyr::filter(tuk, !is.na(p) & p < 0.05) else return(p)
  if (!nrow(tuk)) return(p)
  
  x1 <- match(tuk$group1, lvls); x2 <- match(tuk$group2, lvls)
  span <- abs(x2 - x1)
  p_val <- if ("p.adj" %in% names(tuk)) tuk$p.adj else if ("p" %in% names(tuk)) tuk$p else NA_real_
  ord  <- order(span, dplyr::coalesce(p_val, Inf))
  tuk  <- tuk[ord, , drop = FALSE]; x1 <- x1[ord]; x2 <- x2[ord]
  
  overlaps <- function(a1, a2, b1, b2) !(a2 < b1 || b2 < a1)
  lanes <- list(); lane_id <- integer(nrow(tuk))
  for (i in seq_len(nrow(tuk))) {
    a1 <- min(x1[i], x2[i]); a2 <- max(x1[i], x2[i]); placed <- FALSE
    for (L in seq_along(lanes)) {
      lane_ok <- TRUE
      for (j in seq_len(nrow(lanes[[L]]))) if (overlaps(a1, a2, lanes[[L]]$a1[j], lanes[[L]]$a2[j])) { lane_ok <- FALSE; break }
      if (lane_ok) { lanes[[L]] <- rbind(lanes[[L]], data.frame(a1=a1,a2=a2)); lane_id[i] <- L; placed <- TRUE; break }
    }
    if (!placed) { lanes[[length(lanes)+1]] <- data.frame(a1=a1,a2=a2); lane_id[i] <- length(lanes) }
  }
  
  df_pos <- df %>% dplyr::filter(is.finite(count), count > 0)
  group_max <- tapply(df_pos$count, factor(df_pos$grp, levels = lvls), max, na.rm = TRUE)
  group_max[is.infinite(group_max) | is.na(group_max)] <- 1
  base_h   <- pmax(group_max[tuk$group1], group_max[tuk$group2])
  base_log <- log10(pmax(base_h, 1)); top_log  <- log10(pmax(top_y, 1))
  y_log <- base_log + pmax(pad_dex, min_clr) + (lane_id - 1) * pmax(gap_dex, min_clr)
  cap <- top_log - cap_dex
  
  for (L in unique(lane_id)) {
    idx <- which(lane_id == L); if (length(idx) <= 1) next
    ordL <- idx[order(y_log[idx], na.last = NA)]
    for (k in seq_along(ordL)) {
      i <- ordL[k]
      min_i <- base_log[i] + min_clr + (L - 1) * pmax(gap_dex, min_clr)
      y_log[i] <- max(y_log[i], min_i)
      if (k > 1) y_log[i] <- max(y_log[i], y_log[ordL[k-1]] + min_sep_dex)
    }
    lane_lbl_top <- max(y_log[ordL] + lbl_dex, na.rm = TRUE)
    if (lane_lbl_top > cap) {
      shift <- lane_lbl_top - cap
      y_log[ordL] <- y_log[ordL] - shift
      for (k in seq_along(ordL)) {
        i <- ordL[k]
        min_i <- base_log[i] + min_clr + (L - 1) * pmax(gap_dex, min_clr)
        if (y_log[i] < min_i) y_log[i] <- min_i
        if (k > 1) y_log[i] <- max(y_log[i], y_log[ordL[k-1]] + min_sep_dex)
      }
    }
  }
  
  tip_log <- y_log - tip_dex; lbl_log <- y_log + lbl_dex
  y <- 10^y_log; y0 <- 10^tip_log; ylbl <- 10^lbl_log
  xm <- (pmin(x1, x2) + pmax(x1, x2)) / 2
  lab <- if ("p.adj.signif" %in% names(tuk)) as.character(tuk$p.adj.signif)
  else if ("p.adj" %in% names(tuk)) sprintf("p=%.3g", tuk$p.adj)
  else sprintf("p=%.3g", tuk$p)
  
  seg_df <- data.frame(x1=pmin(x1,x2), x2=pmax(x1,x2), y=y, y0=y0, xm=xm, lab=lab, ylbl=ylbl)
  p +
    geom_segment(data=seg_df, aes(x=x1, xend=x2, y=y, yend=y), linewidth=0.5, inherit.aes=FALSE) +
    geom_segment(data=seg_df, aes(x=x1, xend=x1, y=y0, yend=y), linewidth=0.5, inherit.aes=FALSE) +
    geom_segment(data=seg_df, aes(x=x2, xend=x2, y=y0, yend=y), linewidth=0.5, inherit.aes=FALSE) +
    annotate("text", x=seg_df$xm, y=seg_df$ylbl, label=seg_df$lab, size=3, vjust=0, hjust=0.5)
}

# ---- Crash-guard helpers ----
.as_grob_safe <- function(p) {
  tryCatch(cowplot::as_grob(p),
           error = function(e) cowplot::as_grob(ggplot() + theme_void()))
}
.is_panel_ok <- function(p) inherits(p, c("gg", "grob", "gtable"))

safe_ggsave_pdf <- function(plot, filename, width, height, limitsize = TRUE) {
  tmp <- tempfile(fileext = ".pdf")
  ok <- FALSE
  try({
    ggsave(tmp, plot = plot, width = width, height = height, units = "in",
           dpi = 300, device = grDevices::pdf, useDingbats = FALSE,
           limitsize = limitsize)
    ok <- TRUE
  }, silent = TRUE)
  if (ok) {
    file.copy(tmp, filename, overwrite = TRUE)
  } else {
    message("⚠️ ggsave failed for: ", filename)
  }
  unlink(tmp)
  invisible(ok)
}

# ========= 0) Load DESeq2 objects =========
load("../DESeq2/ofav/host/realModels.RData"); dds_ofav <- dds
load("../DESeq2/ssid/host/realModels.RData"); dds_ssid <- dds
message("ofav colData: ", paste(colnames(SummarizedExperiment::colData(dds_ofav)), collapse=", "))
message("ssid colData: ", paste(colnames(SummarizedExperiment::colData(dds_ssid)), collapse=", "))

# ========= 1) Read lookup and normalize labels =========
df <- read.csv("orthogroup correlation inverse discordant.csv", header = TRUE, check.names = FALSE)
needed <- c("factor","comparison","stress_family","orthogroup_ID",
            "Ofaveolata_ID","Ssiderea_ID","Ofaveolata_lpv","Ssiderea_lpv",
            "relationship","Ofaveolata_annotation","Ssiderea_annotation","gene_name")
stopifnot(all(needed %in% names(df)))

df_norm <- df %>%
  mutate(
    factor = case_when(tolower(factor)=="site" ~ "Site",
                       tolower(factor)=="treatment" ~ "Treatment",
                       TRUE ~ factor),
    relationship = case_when(tolower(relationship)=="direct" ~ "Direct",
                             tolower(relationship)=="inverse" ~ "Inverse",
                             TRUE ~ relationship)
  ) %>%
  mutate(
    Ofaveolata_ID_safe = .safe(Ofaveolata_ID),
    Ssiderea_ID_safe   = .safe(Ssiderea_ID)
  )

# ========= 2) Subsets =========
site_direct       <- df_norm %>% filter(factor=="Site",      relationship=="Direct")  %>% distinct(orthogroup_ID, .keep_all = TRUE)
site_inverse      <- df_norm %>% filter(factor=="Site",      relationship=="Inverse") %>% distinct(orthogroup_ID, .keep_all = TRUE)
treatment_direct  <- df_norm %>% filter(factor=="Treatment", relationship=="Direct")  %>% distinct(orthogroup_ID, .keep_all = TRUE)
treatment_inverse <- df_norm %>% filter(factor=="Treatment", relationship=="Inverse") %>% distinct(orthogroup_ID, .keep_all = TRUE)

cat("\nRows per subset:\n",
    "Site-Direct:   ", nrow(site_direct), "\n",
    "Site-Inverse:  ", nrow(site_inverse), "\n",
    "Treat-Direct:  ", nrow(treatment_direct), "\n",
    "Treat-Inverse: ", nrow(treatment_inverse), "\n", sep = "")

# ========= 3) Export counts CSVs =========
counts_root <- "counts"
dir.create(counts_root, showWarnings = FALSE)
for (s in c(
  "Site_Direct/Ofaveolata","Site_Direct/Ssiderea",
  "Site_Inverse/Ofaveolata","Site_Inverse/Ssiderea",
  "Treatment_Direct/Ofaveolata","Treatment_Direct/Ssiderea",
  "Treatment_Inverse/Ofaveolata","Treatment_Inverse/Ssiderea"
)) dir.create(file.path(counts_root, s), showWarnings = FALSE, recursive = TRUE)

.export_counts_species <- function(dds, gene_id, intgroup, out_dir, species_tag){
  stopifnot(is.character(gene_id), length(gene_id)==1, nzchar(gene_id))
  cd <- SummarizedExperiment::colData(dds)
  if (!(intgroup %in% colnames(cd)))
    stop(sprintf("intgroup '%s' not in colData(%s). Have: %s", intgroup, species_tag, paste(colnames(cd), collapse=", ")))
  if (!(gene_id %in% rownames(dds)))
    stop(sprintf("Gene '%s' not in rownames(%s).", gene_id, species_tag))
  dat <- DESeq2::plotCounts(dds, gene = gene_id, intgroup = intgroup, returnData = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fn <- file.path(out_dir, paste0(.safe(gene_id), ".csv"))
  write.csv(dat, fn, row.names = TRUE)
  message(sprintf("[%s] Wrote: %s (n=%d, cols: %s)", species_tag, normalizePath(fn, mustWork = FALSE), nrow(dat), paste(names(dat), collapse=", ")))
  invisible(fn)
}

export_counts_for_orthogroup <- function(df_sub, og_id, intgroup, out_dir){
  row <- df_sub %>% filter(orthogroup_ID == og_id) %>% slice(1)
  if (nrow(row)==0) stop("OG not found in subset: ", og_id)
  ofav_id <- row$Ofaveolata_ID; ssid_id <- row$Ssiderea_ID
  .export_counts_species(dds_ofav, ofav_id, intgroup, file.path(out_dir, "Ofaveolata"), "Ofaveolata")
  .export_counts_species(dds_ssid, ssid_id, intgroup, file.path(out_dir, "Ssiderea"),  "Ssiderea")
  invisible(NULL)
}

export_counts_batch <- function(df_subset, intgroup, out_dir) {
  if (!nrow(df_subset)) {
    message("Subset is empty; nothing to export.")
    return(invisible(list(success = 0, failed = 0, total = 0)))
  }
  ogs <- unique(df_subset$orthogroup_ID)
  total <- length(ogs)
  message(sprintf("Exporting %d orthogroups to '%s' (intgroup=%s)...",
                  total, normalizePath(out_dir, mustWork = FALSE), intgroup))
  ok <- 0; fail <- 0
  for (og in ogs) {
    tryCatch({
      export_counts_for_orthogroup(df_subset, og, intgroup = intgroup, out_dir = out_dir)
      ok <- ok + 1
    }, error = function(e) {
      message(sprintf("❌ OG %s failed: %s", og, e$message))
      fail <<- fail + 1
    })
  }
  message(sprintf("✅ Export complete for '%s': %d succeeded, %d failed (of %d total).",
                  out_dir, ok, fail, total))
  invisible(list(success = ok, failed = fail, total = total))
}

# Site uses 'site'; Treatment uses 'treat'
res_counts <- list(
  Site_Direct       = export_counts_batch(site_direct,       intgroup = "site",  out_dir = file.path(counts_root,"Site_Direct")),
  Site_Inverse      = export_counts_batch(site_inverse,      intgroup = "site",  out_dir = file.path(counts_root,"Site_Inverse")),
  Treatment_Direct  = export_counts_batch(treatment_direct,  intgroup = "treat", out_dir = file.path(counts_root,"Treatment_Direct")),
  Treatment_Inverse = export_counts_batch(treatment_inverse, intgroup = "treat", out_dir = file.path(counts_root,"Treatment_Inverse"))
)

# ========= 4) Read count CSVs back to long data =========
csv_files <- list.files(counts_root, pattern="\\.csv$", full.names = TRUE, recursive = TRUE)

read_gene_file <- function(file_path){
  dfc <- suppressMessages(readr::read_csv(file_path, show_col_types = FALSE))
  species <- basename(dirname(file_path))
  subset  <- basename(dirname(dirname(file_path)))
  subset_bits <- str_split(subset, "_", simplify = TRUE)
  factor_name <- subset_bits[,1, drop=TRUE]
  relationship <- subset_bits[,2, drop=TRUE]
  gene_id <- tools::file_path_sans_ext(basename(file_path))
  id_col      <- if (species == "Ofaveolata") "Ofaveolata_ID" else "Ssiderea_ID"
  id_col_safe <- paste0(id_col, "_safe")
  
  meta_row <- df_norm %>%
    filter(factor == factor_name,
           relationship == relationship,
           (!!rlang::sym(id_col) == gene_id) | (!!rlang::sym(id_col_safe) == gene_id)) %>% slice(1)
  if (!nrow(meta_row)) {
    meta_row <- df_norm %>%
      filter((!!rlang::sym(id_col) == gene_id) | (!!rlang::sym(id_col_safe) == gene_id)) %>% slice(1)
  }
  
  tibble(
    species       = species,
    gene_id       = gene_id,
    orthogroup_ID = meta_row$orthogroup_ID %||% NA_character_,
    gene_name     = meta_row$gene_name %||% NA_character_,
    stress_family = meta_row$stress_family %||% NA_character_,
    relationship  = relationship,
    factor        = factor_name,
    subset        = paste(factor_name, relationship, sep = "_"),
    count_data    = list(dfc)
  )
}

counts_list <- purrr::map_dfr(csv_files, read_gene_file)
counts_long <- counts_list %>% tidyr::unnest(count_data)

# ========= 5) Plotting + stats =========
species_titles <- list(
  Ofaveolata = expression(italic("O. faveolata")),
  Ssiderea   = expression(italic("S. siderea"))
)
site_levels   <- c("Emerald","Rainbow","Star","MacN")
site_colors   <- c("#018571","#80cdc1","#dfc27d","#a6611a")
treat_levels  <- c("CC","LC","CH","LH")
treat_colors  <- c("#92c5de","#0571b0","#f4a582","#ca0020")
stress_family_palette <- c(
  "Apoptosis/cell death"          = "#EF4444",
  "DNA repair"                    = "#10B981",
  "Heat shock proteins"           = "#8B5CF6",
  "Immune response"               = "#F59E0B",
  "Oxidative stress/antioxidants" = "#3B82F6",
  "Ubiquitin-proteasome"          = "#6366F1",
  "Other stress-related"          = "#6B7280"
)
use_fixed_y_limits <- TRUE
fixed_y_limits <- c(3, 12000)

site_labels <- c(
  "Emerald" = "Emerald\nReef",
  "Rainbow" = "Rainbow\nReef",
  "Star"    = "Star\nIsland",
  "MacN"    = "MacArthur\nNorth"
)
treat_labels <- c(
  "CC" = "Contemporary\npH + Ambient\nTemperature",
  "LC" = "Acidified\n+ Ambient\nTemperature",
  "CH" = "Contemporary\npH +\nBleaching",
  "LH" = "Acidified\n+\nBleaching"
)

panel_with_colored_title <- function(left_plot, right_plot, title_text, stress_family){
  fill_col <- if (!is.null(stress_family) && stress_family %in% names(stress_family_palette)) {
    stress_family_palette[stress_family]
  } else "#6B7280"
  
  gL <- ggplotGrob(left_plot); gR <- ggplotGrob(right_plot)
  maxw <- grid::unit.pmax(gL$widths[2:5], gR$widths[2:5])
  gL$widths[2:5] <- maxw; gR$widths[2:5] <- maxw
  
  core <- cowplot::plot_grid(cowplot::as_grob(gL), cowplot::as_grob(gR), ncol = 2, rel_widths = c(1, 0.92))
  
  title_strip <- ggplot() +
    geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf), fill=fill_col, color=NA) +
    annotate("text", x=0, y=0, label=title_text, color="white", fontface="bold", size=4.2, hjust=0.5, vjust=0.5) +
    coord_cartesian(xlim=c(-1,1), ylim=c(-1,1), expand=FALSE) + theme_void()
  
  cowplot::plot_grid(title_strip, core, ncol=1, rel_heights=c(0.12, 1))
}

compute_stats <- function(df_grp){
  df2 <- df_grp %>% dplyr::filter(!is.na(grp), is.finite(count))
  out <- list(aov_text = NULL, aov_table = NULL, tukey = NULL, aov_label = NULL, aov_p = NA_real_)
  if (!nrow(df2)) return(out)
  
  lvl_counts <- table(droplevels(df2$grp))
  if (sum(lvl_counts > 0) < 2L) return(out)
  
  aov_res <- tryCatch(aov(count ~ grp, data = df2), error = function(e) NULL)
  got_aov <- !is.null(aov_res)
  pv <- NA_real_; Fv <- NA_real_; df1 <- NA_integer_; df2res <- NA_integer_
  
  if (got_aov) {
    sm <- summary(aov_res)[[1]]
    out$aov_table <- as.data.frame(sm) %>% tibble::rownames_to_column("Term") %>%
      dplyr::mutate(across(where(is.numeric), ~ round(.x, 4)))
    term_names <- attr(terms(aov_res), "term.labels")
    term <- if (length(term_names)) term_names[1] else "grp"
    row_idx <- which(rownames(sm) == term)
    if (!length(row_idx)) row_idx <- which(grepl(paste0("^", term, "$"), rownames(sm)))
    if (length(row_idx)) {
      df1     <- sm[row_idx, "Df"]
      df2res  <- sm["Residuals", "Df"]
      Fv      <- sm[row_idx, "F value"]
      pv      <- sm[row_idx, "Pr(>F)"]
    }
  }
  
  if (!is.finite(pv)) {
    welch <- tryCatch(oneway.test(count ~ grp, data = df2, var.equal = FALSE), error = function(e) NULL)
    if (!is.null(welch)) {
      pv <- as.numeric(welch$p.value)
      df1 <- unname(signif(welch$parameter[1], 3))
      df2res <- unname(signif(welch$parameter[2], 3))
      Fv <- unname(signif(welch$statistic[[1]], 3))
    }
  }
  if (!is.finite(pv)) {
    kw <- tryCatch(kruskal.test(count ~ grp, data = df2), error = function(e) NULL)
    if (!is.null(kw)) {
      pv <- as.numeric(kw$p.value)
      Fv <- NA_real_; df1 <- NA; df2res <- NA
    }
  }
  
  out$aov_p <- pv
  
  # --- Plain text label (CSV/logs) with star when p<0.05 ---
  pv_str_plain <- if (is.finite(pv)) formatC(pv, format = "g", digits = 3) else "NA"
  star_plain   <- if (is.finite(pv) && pv < 0.05) "*" else ""
  if (is.finite(Fv) && is.finite(pv) && is.finite(as.numeric(df1)) && is.finite(as.numeric(df2res))) {
    out$aov_text <- sprintf("ANOVA: F(%s,%s)=%.2f, p=%s%s", df1, df2res, as.numeric(Fv), pv_str_plain, star_plain)
  } else if (is.finite(pv)) {
    out$aov_text <- sprintf("ANOVA: p=%s%s", pv_str_plain, star_plain)
  } else {
    out$aov_text <- "ANOVA: p=NA"
  }
  
  # --- Plotmath label with star when p<0.05 ---
  if (is.finite(pv)) {
    pretty_p <- if (pv < 1e-4) "'< 1e-4'"
    else if (pv < 0.001) "'< 0.001'"
    else formatC(pv, format = "g", digits = 3)
    star_plot <- if (pv < 0.05) ", '*'" else ""  # appended as another paste() arg
    
    Fv_str <- if (is.finite(Fv)) formatC(as.numeric(Fv), format = "f", digits = 2) else "NA"
    if (is.finite(Fv) && is.finite(as.numeric(df1)) && is.finite(as.numeric(df2res))) {
      out$aov_label <- sprintf(
        "paste(italic(F)[%s*','*%s], ' = ', %s, ', ', italic(p), ' = ', %s%s)",
        df1, df2res, Fv_str, pretty_p, star_plot
      )
    } else {
      out$aov_label <- sprintf("paste(italic(p), ' = ', %s%s)", pretty_p, star_plot)
    }
  } else {
    out$aov_label <- NULL
  }
  
  # Tukey only makes sense for classic AOV; compute but we'll gate display elsewhere
  out$tukey <- if (got_aov) tryCatch(rstatix::tukey_hsd(aov_res), error = function(e) NULL) else NULL
  out
}

pair_order <- function(levels_vec){
  pairs <- list()
  for (i in seq_along(levels_vec)) for (j in seq(i+1, length(levels_vec))) pairs[[length(pairs)+1]] <- c(levels_vec[i], levels_vec[j])
  do.call(rbind, pairs) %>% as.data.frame() %>% setNames(c("group1","group2"))
}
order_tukey_by_levels <- function(tukey_tbl, levels_vec){
  if (is.null(tukey_tbl) || !nrow(tukey_tbl)) return(NULL)
  po <- pair_order(levels_vec)
  tukey_tbl %>%
    mutate(order_key = match(paste(group1, group2, sep="__"), paste(po$group1, po$group2, sep="__"))) %>%
    arrange(order_key, p.adj) %>% select(-order_key)
}

make_boxplot <- function(df, factor_name, species_label,
                         species_name,
                         tukey_tbl = NULL, aov_label = NULL, aov_text = NULL,
                         y_limits = NULL,
                         suppress_y_title = FALSE,   # NEW
                         suppress_x_text  = FALSE) { # NEW
  # determine levels, colors, and labels
  if (factor_name == "Site") {
    lvls <- site_levels
    cols <- site_colors
    label_map <- site_labels
  } else {
    lvls <- treat_levels
    cols <- treat_colors
    label_map <- treat_labels
  }
  df$grp <- factor(df$grp, levels = lvls)
  
  # base plot
  p <- ggpubr::ggboxplot(df, x = "grp", y = "count", color = "grey30", fill = "grp",
                         add = "jitter", add.params = list(size = 1, jitter = 0.25),
                         width = 0.7, size = 0.5) +
    labs(
      title = species_label,
      x = NULL,
      y = if (species_name == "Ssiderea") NULL else "Log Normalized Counts",
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(
        size = if (factor_name == "Treatment") 8 else 10,
        lineheight = 0.9,
        angle = 0, hjust = 0.5, vjust = 1
      ),
      axis.text.y  = if (species_name == "Ssiderea") element_blank() else element_text(),
      axis.ticks.y = if (species_name == "Ssiderea") element_blank() else element_line(),
      plot.margin = margin(t = 6, r = 6, b = 18, l = 6),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      panel.grid = element_blank()
    ) +
    scale_fill_manual(values = setNames(cols, lvls)) +
    scale_x_discrete(labels = label_map)
  
  # --- unified y scale ---
  if (!is.null(y_limits)) {
    ylims_used <- y_limits
  } else {
    pos <- df$count[is.finite(df$count) & df$count > 0]
    if (length(pos)) {
      min_pos <- min(pos, na.rm = TRUE)
      max_pos <- max(pos, na.rm = TRUE)
      lower_pad <- 1/1.3
      upper_pad <- 1.3
      ymin <- max(min_pos * lower_pad, min_pos * 0.5, .Machine$double.eps)
      ymax <- max_pos * upper_pad
    } else {
      ymin <- 0.1; ymax <- 1
    }
    ylims_used <- c(ymin, ymax)
  }
  
  # ensure room for Tukey if present
  if (!is.null(tukey_tbl) && nrow(tukey_tbl)) {
    need_top <- .bracket_required_top(
      df, tukey_tbl, lvls,
      pad_dex   = 0.15,
      gap_dex   = 0.15,
      lbl_dex   = 0.005,
      min_sep_dex = 0.08,
      micro_dex = 0.020
    )
    if (is.finite(need_top)) ylims_used[2] <- max(ylims_used[2], need_top)
  }
  
  p <- p +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = expansion(mult = c(0.02, 0.12))
    ) +
    coord_cartesian(ylim = ylims_used, clip = "off")
  
  # --- Tukey brackets ---
  top_y <- ylims_used[2] * 0.98
  if (!is.null(tukey_tbl) && nrow(tukey_tbl)) {
    p <- .add_sig_brackets(
      p, tukey_tbl, lvls, top_y, df,
      pad_dex = 0.2,
      gap_dex = 0.15,
      tip_dex = 0.025,
      lbl_dex = 0.005
    )
  }
  
  # --- ANOVA label (bottom-center) ---
  lab_text  <- if (!is.null(aov_label) && nzchar(aov_label)) aov_label else aov_text
  use_parse <- !is.null(aov_label) && nzchar(aov_label)
  center_x  <- mean(seq_along(lvls))
  bottom_y  <- 10^(log10(ylims_used[1]) - 0.03)
  if (!is.null(lab_text) && nzchar(lab_text)) {
    p <- p + annotate("text", x = center_x, y = bottom_y, label = lab_text,
                      size = 3.3, vjust = 0, hjust = 0.5, parse = use_parse)
  }
  
  # --- NEW: per-position axis suppression for one-pagers ---
  if (isTRUE(suppress_y_title)) p <- p + theme(axis.title.y = element_blank())
  if (isTRUE(suppress_x_text))  p <- p + theme(axis.text.x  = element_blank())
  
  list(plot = p)
}

# ========= 5b) One-page composite helpers (8.5x11 auto layout) =========

# (keep .best_grid_dims and .page_title_strip you already have, or paste from below if missing)

.best_grid_dims <- function(n, content_w, content_h, target_ar, max_rc = 20L) {
  cand <- expand.grid(r = 1:max_rc, c = 1:max_rc)
  cand <- cand[cand$r * cand$c >= n, , drop = FALSE]
  cand$cell_w <- content_w / cand$c
  cand$cell_h <- content_h / cand$r
  cand$ar     <- cand$cell_w / cand$cell_h
  cand$score  <- abs(cand$ar - target_ar)
  cand$empty  <- (cand$r * cand$c) - n
  cand$area   <- cand$cell_w * cand$cell_h
  cand <- cand[order(cand$score, cand$empty, -cand$area), ]
  head(cand, 1)
}

.page_title_strip <- function(title_text) {
  ggplot() +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), fill="grey10", color=NA) +
    annotate("text", x=0, y=0, label=title_text, color="white", fontface="bold", size=12, hjust=0.5, vjust=0.5) +
    coord_cartesian(xlim=c(-1,1), ylim=c(-1,1), expand=FALSE) +
    theme_void()
}

.build_stress_family_legend <- function(palette, ncol = 3, title = "Stress family") {
  df_leg <- data.frame(fam = factor(names(palette), levels = names(palette)))
  p_leg <- ggplot(df_leg, aes(x = fam, y = 1, fill = fam)) +
    geom_col() +
    scale_fill_manual(values = palette, name = title) +
    guides(fill = guide_legend(ncol = ncol, byrow = TRUE)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text  = element_text(size = 8),
      legend.key.height = unit(0.4, "lines"),
      legend.key.width  = unit(0.8, "lines")
    )
  cowplot::get_legend(p_leg)
}

# ===== Fixed-grid, no-shrink, landscape, legend-in-empty-cell exporter =====

# ===== Fixed-grid, no-shrink, landscape, legend-in-empty-cell (NO page title) =====
export_subset_fixed_grid_no_shrink <- function(dat_sub, fac, rel, ogs, out_dir,
                                               rows, cols,
                                               panel_w = 8, panel_h = 4.5,
                                               margin_l = 0.35, margin_r = 0.35,
                                               margin_t = 0.30, margin_b = 0.30,
                                               legend_ncol = 3) {
  if (!length(ogs)) return(invisible(NULL))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Bigger legend
  .build_stress_family_legend_big <- function(palette, ncol = 3, title = "Stress family") {
    df_leg <- data.frame(fam = factor(names(palette), levels = names(palette)))
    p_leg <- ggplot(df_leg, aes(x = fam, y = 1, fill = fam)) +
      geom_col() +
      scale_fill_manual(values = palette, name = title) +
      guides(fill = guide_legend(ncol = ncol, byrow = TRUE)) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 11),
        legend.key.height = unit(0.6, "lines"),
        legend.key.width  = unit(1.2, "lines"),
        legend.box.margin = margin(2,2,2,2)
      )
    cowplot::get_legend(p_leg)
  }
  
  final_panels <- list()
  lvls <- if (fac == "Site") site_levels else treat_levels
  
  max_slots <- rows * cols
  n_place <- min(length(ogs), max_slots)
  
  # First pass: collect per-OG data
  tmp <- vector("list", n_place)
  for (i in seq_len(n_place)) {
    og <- ogs[i]
    row_meta <- dat_sub %>% dplyr::filter(orthogroup_ID==og) %>% dplyr::slice(1)
    if (!nrow(row_meta)) next
    of_df <- dat_sub %>% dplyr::filter(orthogroup_ID==og, species=="Ofaveolata")
    ss_df <- dat_sub %>% dplyr::filter(orthogroup_ID==og, species=="Ssiderea")
    if (!nrow(of_df) || !nrow(ss_df)) next
    pos <- c(of_df$count, ss_df$count); pos <- pos[is.finite(pos) & pos > 0]
    if (length(pos)) { ymin <- max(min(pos)/1.3, min(pos)*0.5, .Machine$double.eps); ymax <- max(pos)*1.3 } else { ymin <- 0.1; ymax <- 1 }
    ylims <- c(ymin, ymax)
    of_stats <- compute_stats(of_df)
    ss_stats <- compute_stats(ss_df)
    tmp[[i]] <- list(of_df=of_df, ss_df=ss_df, ylims=ylims, of_stats=of_stats, ss_stats=ss_stats,
                     gname = row_meta$gene_name %>% .[1], sfam = row_meta$stress_family %>% .[1])
  }
  
  # Second pass: build grobs with correct axis suppression
  for (i in seq_len(n_place)) {
    it <- tmp[[i]]; if (is.null(it)) next
    # Grid position
    r <- ceiling(i / cols)
    c <- ((i - 1) %% cols) + 1
    # X labels on any plot that does NOT have another real plot under it
    has_below <- (i + cols) <= n_place
    suppress_x_text <- has_below
    # Y title only on first column
    suppress_y_title_left <- !(c == 1)
    
    of_tukey_tbl <- if (!is.na(it$of_stats$aov_p) && is.finite(it$of_stats$aov_p) && it$of_stats$aov_p < 0.05)
      order_tukey_by_levels(it$of_stats$tukey, lvls) else NULL
    ss_tukey_tbl <- if (!is.na(it$ss_stats$aov_p) && is.finite(it$ss_stats$aov_p) && it$ss_stats$aov_p < 0.05)
      order_tukey_by_levels(it$ss_stats$tukey, lvls) else NULL
    
    of_plot <- make_boxplot(it$of_df, fac, species_titles$Ofaveolata,
                            species_name = "Ofaveolata",
                            tukey_tbl = of_tukey_tbl,
                            aov_label = it$of_stats$aov_label,
                            aov_text  = it$of_stats$aov_text,
                            y_limits  = it$ylims,
                            suppress_y_title = suppress_y_title_left,
                            suppress_x_text  = suppress_x_text)$plot
    
    ss_plot <- make_boxplot(it$ss_df, fac, species_titles$Ssiderea,
                            species_name = "Ssiderea",
                            tukey_tbl = ss_tukey_tbl,
                            aov_label = it$ss_stats$aov_label,
                            aov_text  = it$ss_stats$aov_text,
                            y_limits  = it$ylims,
                            suppress_y_title = TRUE,
                            suppress_x_text  = suppress_x_text)$plot
    
    # Keep your colored stress-family bar inside each mini-panel:
    final_panels[[length(final_panels)+1L]] <- panel_with_colored_title(of_plot, ss_plot,
                                                                        paste0(it$gname,"  (",ogs[i],")"),
                                                                        it$sfam)
  }
  
  # Fill remaining slots; put BIG legend in the first empty slot
  total_now <- length(final_panels)
  empty <- max_slots - total_now
  if (empty > 0) {
    legend_panel <- {
      leg <- .build_stress_family_legend_big(stress_family_palette, ncol = legend_ncol)
      cowplot::ggdraw() + cowplot::draw_grob(leg, 0.5, 0.5, 0.98, 0.98)
    }
    final_panels[[length(final_panels)+1L]] <- legend_panel
    if (empty > 1) {
      blanks <- replicate(empty - 1L, ggplot() + theme_void(), simplify = FALSE)
      final_panels <- c(final_panels, blanks)
    }
  }
  
  # Grid only (NO page title) — sanitize panels first
  final_panels <- Filter(.is_panel_ok, final_panels)
  final_panels <- lapply(final_panels, .as_grob_safe)
  
  grid_core <- cowplot::plot_grid(
    plotlist = final_panels,
    nrow = rows, ncol = cols,
    align = "hv", axis = "tblr"
  )
  
  # Page size (landscape), no title height
  page_w <- margin_l + cols * panel_w + margin_r
  page_h <- margin_t + rows * panel_h + margin_b
  
  outfile <- file.path(out_dir, sprintf("ONEPAGE_%s_%s_fixedgrid.pdf", fac, rel))
  invisible(safe_ggsave_pdf(grid_core, outfile, width = page_w, height = page_h, limitsize = FALSE))
  
  # Aggressive cleanup to avoid OOM before next subset
  rm(final_panels, grid_core); gc()
  invisible(outfile)
}

# ========= Output root =========
panels_root <- "plots_panels"
dir.create(panels_root, showWarnings = FALSE)

targets <- tibble::tribble(
  ~factor_name, ~relationship, ~out_subdir,
  "Site", "Direct", "Site_Direct",
  "Site", "Inverse","Site_Inverse",
  "Treatment","Direct","Treatment_Direct",
  "Treatment","Inverse","Treatment_Inverse"
)

# ========= Main loop (with fixed-grid, no-shrink, legend-in-empty-cell) =========
for (k in seq_len(nrow(targets))) {
  fac <- targets$factor_name[k]
  rel <- targets$relationship[k]
  out_dir <- file.path(panels_root, targets$out_subdir[k])
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "tukey"), showWarnings = FALSE)
  dir.create(file.path(out_dir, "anova"), showWarnings = FALSE)
  
  dat_sub <- counts_long %>% dplyr::filter(factor == fac, relationship == rel)
  if (!nrow(dat_sub)) {
    message(sprintf("No rows for %s | %s", fac, rel))
    next
  }
  
  if (fac == "Site") {
    stopifnot("site" %in% names(dat_sub))
    dat_sub <- dat_sub %>% dplyr::mutate(grp = factor(site, levels = site_levels))
    lvls <- site_levels
  } else {
    stopifnot("treat" %in% names(dat_sub))
    dat_sub <- dat_sub %>% dplyr::mutate(grp = factor(treat, levels = treat_levels))
    lvls <- treat_levels
  }
  
  # Order by stress_family then orthogroup_ID for panel sequence
  ogs <- df_norm %>%
    dplyr::filter(factor==fac, relationship==rel) %>%
    dplyr::distinct(orthogroup_ID, .keep_all = TRUE) %>%
    dplyr::arrange(stress_family, orthogroup_ID) %>%
    dplyr::pull(orthogroup_ID)
  ogs <- intersect(ogs, unique(dat_sub$orthogroup_ID))
  
  panel_list <- list()
  p_i <- 0L
  
  for (og in ogs) {
    row_meta <- dat_sub %>% dplyr::filter(orthogroup_ID == og) %>% dplyr::slice(1)
    gname <- row_meta$gene_name %>% .[1]
    sfam  <- row_meta$stress_family %>% .[1]
    
    of_df <- dat_sub %>% dplyr::filter(orthogroup_ID == og, species == "Ofaveolata")
    ss_df <- dat_sub %>% dplyr::filter(orthogroup_ID == og, species == "Ssiderea")
    
    if (!nrow(of_df) || !nrow(ss_df)) {
      missing_species <- paste(
        c("Ofaveolata", "Ssiderea")[c(!nrow(of_df), !nrow(ss_df))],
        collapse = ", "
      )
      message(sprintf("Skipping OG=%s (missing species: %s)", og, missing_species))
      next
    }
    
    # shared y-limits across both species for this OG (log-safe)
    counts_both <- c(of_df$count, ss_df$count)
    pos <- counts_both[is.finite(counts_both) & counts_both > 0]
    
    if (length(pos)) {
      min_pos <- min(pos, na.rm = TRUE)
      max_pos <- max(pos, na.rm = TRUE)
      lower_pad <- 1 / 1.3
      upper_pad <- 1.3
      ymin <- max(min_pos * lower_pad, min_pos * 0.5, .Machine$double.eps)
      ymax <- max_pos * upper_pad
    } else {
      ymin <- 0.1
      ymax <- 1
    }
    ylims <- c(ymin, ymax)
    
    of_stats <- compute_stats(of_df)
    ss_stats <- compute_stats(ss_df)
    
    # Save stats CSVs
    if (!is.null(of_stats$tukey) && nrow(of_stats$tukey))
      readr::write_csv(order_tukey_by_levels(of_stats$tukey, lvls),
                       file.path(out_dir, "tukey", paste0("Tukey_Ofaveolata_", .safe(og), ".csv")))
    if (!is.null(ss_stats$tukey) && nrow(ss_stats$tukey))
      readr::write_csv(order_tukey_by_levels(ss_stats$tukey, lvls),
                       file.path(out_dir, "tukey", paste0("Tukey_Ssiderea_", .safe(og), ".csv")))
    if (!is.null(of_stats$aov_table))
      readr::write_csv(of_stats$aov_table,
                       file.path(out_dir, "anova", paste0("ANOVA_Ofaveolata_", .safe(og), ".csv")))
    if (!is.null(ss_stats$aov_table))
      readr::write_csv(ss_stats$aov_table,
                       file.path(out_dir, "anova", paste0("ANOVA_Ssiderea_", .safe(og), ".csv")))
    
    # Gate Tukey: only pass it through if omnibus ANOVA is significant
    of_tukey_tbl <- if (!is.na(of_stats$aov_p) && is.finite(of_stats$aov_p) && of_stats$aov_p < 0.05)
      order_tukey_by_levels(of_stats$tukey, lvls) else NULL
    ss_tukey_tbl <- if (!is.na(ss_stats$aov_p) && is.finite(ss_stats$aov_p) && ss_stats$aov_p < 0.05)
      order_tukey_by_levels(ss_stats$tukey, lvls) else NULL
    
    of_plot <- make_boxplot(of_df, fac, species_titles$Ofaveolata,
                            species_name = "Ofaveolata",
                            tukey_tbl = of_tukey_tbl,
                            aov_label = of_stats$aov_label,
                            aov_text  = of_stats$aov_text,
                            y_limits  = ylims)$plot
    
    ss_plot <- make_boxplot(ss_df, fac, species_titles$Ssiderea,
                            species_name = "Ssiderea",
                            tukey_tbl = ss_tukey_tbl,
                            aov_label = ss_stats$aov_label,
                            aov_text  = ss_stats$aov_text,
                            y_limits  = ylims)$plot
    
    panel <- panel_with_colored_title(of_plot, ss_plot, paste0(gname, "  (", og, ")"), sfam)
    p_i <- p_i + 1L
    outfile <- file.path(out_dir,
                         paste0("PANEL_", sprintf("%03d", p_i), "_", .safe(og), "_", .safe(gname), ".pdf"))
    # old:
    # ggsave(outfile, panel, width = 8, height = 4.5, dpi = 300, device = grDevices::pdf, useDingbats = FALSE)
    # new:
    safe_ggsave_pdf(panel, outfile, width = 8, height = 4.5, limitsize = TRUE)
    message("Wrote panel: ", normalizePath(outfile, mustWork = FALSE))
    panel_list[[p_i]] <- panel
  }
  
  if (length(panel_list)) {
    # (Optional) keep your ALL_panels export
    ggexport(plotlist = panel_list,
             filename = file.path(out_dir, paste0("ALL_panels_", fac, "_", rel, ".pdf")),
             nrow = 1, ncol = 1, width = 8, height = 4.5, dpi = 300)
    
    # Fixed-grid, no-shrink, legend-in-empty-cell (LANDSCAPE), NO page title
    if (fac == "Site" && rel == "Direct") {
      export_subset_fixed_grid_no_shrink(dat_sub, fac, rel, ogs, out_dir, rows = 5, cols = 4)
    } else if (fac == "Site" && rel == "Inverse") {
      export_subset_fixed_grid_no_shrink(dat_sub, fac, rel, ogs, out_dir, rows = 4, cols = 4)  # your choice
    } else if (fac == "Treatment" && rel == "Direct") {
      export_subset_fixed_grid_no_shrink(dat_sub, fac, rel, ogs, out_dir, rows = 8, cols = 7)
    } else if (fac == "Treatment" && rel == "Inverse") {
      export_subset_fixed_grid_no_shrink(dat_sub, fac, rel, ogs, out_dir, rows = 6, cols = 5)
    }
  }
  
  # ✅ per-subset panel summary
  message(sprintf("📊 %d panels exported for %s | %s → %s",
                  length(panel_list), fac, rel, normalizePath(out_dir, mustWork = FALSE)))
  
  # Free memory before moving to the next subset
  rm(dat_sub, panel_list, of_plot, ss_plot, of_df, ss_df, of_stats, ss_stats, ylims, counts_both, pos)
  gc()
  
}

# ========= Summaries =========
message("\nSummary by subset (panels):")
for (k in seq_len(nrow(targets))) {
  subdir <- file.path("plots_panels", targets$out_subdir[k])
  npan <- length(list.files(subdir, pattern="^PANEL_.*\\.pdf$", full.names=TRUE))
  message(sprintf("  %s | %s : %d panels", targets$factor_name[k], targets$relationship[k], npan))
}

message("\nCounts export summary (from earlier):"); print(res_counts)
message("\n✅ Done. Counts in 'counts/', panels + Tukey/ANOVA CSVs in 'plots_panels/'.")

message("\nSummary by subset (Tukey & ANOVA CSVs):")
for (k in seq_len(nrow(targets))) {
  subdir  <- file.path("plots_panels", targets$out_subdir[k])
  tuk_dir <- file.path(subdir, "tukey")
  aov_dir <- file.path(subdir, "anova")
  ntukey <- if (dir.exists(tuk_dir)) length(list.files(tuk_dir, pattern="\\.csv$", full.names=TRUE)) else 0L
  nanova <- if (dir.exists(aov_dir)) length(list.files(aov_dir, pattern="\\.csv$", full.names=TRUE)) else 0L
  message(sprintf("  %s | %s : Tukey=%d, ANOVA=%d", targets$factor_name[k], targets$relationship[k], ntukey, nanova))
}

# 🔎 List the single-page outputs
onepagers <- list.files("plots_panels", pattern = "^ONEPAGE_.*\\.pdf$", recursive = TRUE, full.names = TRUE)
message("\nFound ", length(onepagers), " one-page PDFs.")
if (length(onepagers)) cat(paste0(" - ", normalizePath(onepagers, mustWork = FALSE), "\n"), sep = "")
