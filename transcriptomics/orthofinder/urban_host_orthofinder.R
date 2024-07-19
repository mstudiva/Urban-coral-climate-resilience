#### PACKAGES ####

library(tidyverse)
library(VennDiagram)
library(pheatmap)
library(reshape2)
library(RColorBrewer)


#### ORTHOFINDER ####

# Install orthofinder on your local machine using the tutorials (https://davidemms.github.io/menu/tutorials.html)

# Copy your translated protein fasta files (_out_PRO.fas) that you want to compare into a directory called 'orthofinder'
# If you have not already filtered by the longest contig per isogroup (by using fasta2SBH.pl during transcriptome annotation), follow step 7 of tutorial 2 above

# Run the following command in Terminal: 'orthofinder -f orthofinder/'
# Check the number of genes assigned to orthogroups (e.g., 'OrthoFinder assigned 48479 genes (83.8% of total) to 13598 orthogroups')
# Ideally, it should be >80%


#### ORTHOLOGS ####

load("orthofinder_DEGs.RData") # if previously run
orthologs <- read.table(file = "orthofinder/OrthoFinder/Results_Oct08/Orthologues/Orthologues_Ofaveolata_out_PRO/Ofaveolata_out_PRO__v__Siderastrea_out_PRO.tsv", sep = "\t", header = TRUE, quote="", fill=FALSE)


#### DESEQ IMPORT TREATMENT ####

ofav_LC_CC_lpv <- read.csv(file = "../../DESeq2/ofav/host/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_CH_CC_lpv <- read.csv(file = "../../DESeq2/ofav/host/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_CC_lpv <- read.csv(file = "../../DESeq2/ofav/host/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_CH_LC_lpv <- read.csv(file = "../../DESeq2/ofav/host/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_CH_lpv <- read.csv(file = "../../DESeq2/ofav/host/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_LH_LC_lpv <- read.csv(file = "../../DESeq2/ofav/host/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ssid_LC_CC_lpv <- read.csv(file = "../../DESeq2/ssid/host/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_CH_CC_lpv <- read.csv(file = "../../DESeq2/ssid/host/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_CC_lpv <- read.csv(file = "../../DESeq2/ssid/host/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_CH_LC_lpv <- read.csv(file = "../../DESeq2/ssid/host/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_CH_lpv <- read.csv(file = "../../DESeq2/ssid/host/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_LH_LC_lpv <- read.csv(file = "../../DESeq2/ssid/host/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)


#### DESEQ IMPORT SITE ####
ofav_Rainbow_Emerald_lpv <- read.csv(file = "../../DESeq2/ofav/host/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_Star_Emerald_lpv <- read.csv(file = "../../DESeq2/ofav/host/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Emerald_lpv <- read.csv(file = "../../DESeq2/ofav/host/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_Star_Rainbow_lpv <- read.csv(file = "../../DESeq2/ofav/host/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Rainbow_lpv <- read.csv(file = "../../DESeq2/ofav/host/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ofav_MacN_Star_lpv <- read.csv(file = "../../DESeq2/ofav/host/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ofav" = lpv)

ssid_Rainbow_Emerald_lpv <- read.csv(file = "../../DESeq2/ssid/host/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_Star_Emerald_lpv <- read.csv(file = "../../DESeq2/ssid/host/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Emerald_lpv <- read.csv(file = "../../DESeq2/ssid/host/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_Star_Rainbow_lpv <- read.csv(file = "../../DESeq2/ssid/host/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Rainbow_lpv <- read.csv(file = "../../DESeq2/ssid/host/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)

ssid_MacN_Star_lpv <- read.csv(file = "../../DESeq2/ssid/host/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_ssid" = lpv)


#### DEG MATCHING TREATMENT ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds ofav and ssid gene annotations, and 5) then pulls on corresponding KOG classes

# LC vs CC for both species
orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LC_CC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LC_CC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_CH_CC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_CH_CC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LH_CC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LH_CC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_CH_LC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_CH_LC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LH_CH_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LH_CH_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_LH_LC_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_LH_LC_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_Rainbow_Emerald_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_Rainbow_Emerald_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_Star_Emerald_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_Star_Emerald_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_MacN_Emerald_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_MacN_Emerald_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_Star_Rainbow_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_Star_Rainbow_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_MacN_Rainbow_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_MacN_Rainbow_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  inner_join(ssid_MacN_Star_lpv, by = c("Protein_ssid" = "gene")) %>%
  inner_join(ofav_MacN_Star_lpv, by = c("Protein_ofav" = "gene")) %>%
  filter(abs(lpv_ssid) >= 1 & abs(lpv_ofav) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_ssid = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ssid" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ofav/young/Ofaveolata_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_ofav = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_ofav" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/ssid/magana/Siderastrea_iso2kogClass.tab",
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


#### VENN DIAGRAMS TREATMENT ####

# first creating a set of up/downregulated DEGs by species
LC_CC %>%
  filter(lpv_ssid >= 1) %>%
  pull(Orthogroup) -> ssid_up

LC_CC %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

LC_CC %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

LC_CC %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
pdf(file="Venn_LC_CC.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
CH_CC %>%
  filter(lpv_ssid >= 1) %>%
  pull(Orthogroup) -> ssid_up

CH_CC %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

CH_CC %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

CH_CC %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
  pull(Orthogroup) -> ssid_up

LH_CC %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

LH_CC %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

LH_CC %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
  pull(Orthogroup) -> ssid_up

CH_LC %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

CH_LC %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

CH_LC %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
LH_CH %>%
  filter(lpv_ssid >= 1) %>%
  pull(Orthogroup) -> ssid_up

LH_CH %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

LH_CH %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

LH_CH %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
pdf(file="Venn_LH_CH.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
LH_LC %>%
  filter(lpv_ssid >= 1) %>%
  pull(Orthogroup) -> ssid_up

LH_LC %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

LH_LC %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

LH_LC %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
  pull(Orthogroup) -> ssid_up

Rainbow_Emerald %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

Rainbow_Emerald %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

Rainbow_Emerald %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
  pull(Orthogroup) -> ssid_up

Star_Emerald %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

Star_Emerald %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

Star_Emerald %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
  pull(Orthogroup) -> ssid_up

MacN_Emerald %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

MacN_Emerald %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

MacN_Emerald %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
  pull(Orthogroup) -> ssid_up

Star_Rainbow %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

Star_Rainbow %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

Star_Rainbow %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
MacN_Star %>%
  filter(lpv_ssid >= 1) %>%
  pull(Orthogroup) -> ssid_up

MacN_Star %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

MacN_Star %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

MacN_Star %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
pdf(file="Venn_MacN_Star.pdf", height=10, width=12)
grid.draw(venn)
dev.off()

# first creating a set of up/downregulated DEGs by species
MacN_Rainbow %>%
  filter(lpv_ssid >= 1) %>%
  pull(Orthogroup) -> ssid_up

MacN_Rainbow %>%
  filter(lpv_ssid <= -1) %>%
  pull(Orthogroup) -> ssid_down

MacN_Rainbow %>%
  filter(lpv_ofav >= 1) %>%
  pull(Orthogroup) -> ofav_up

MacN_Rainbow %>%
  filter(lpv_ofav <= -1) %>%
  pull(Orthogroup) -> ofav_down

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
LC_CC %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  LC_CC_heatmap

CH_CC %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  CH_CC_heatmap

LH_CC %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  LH_CC_heatmap

CH_LC %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  CH_LC_heatmap

LH_CH %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  LH_CH_heatmap

LH_LC %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  LH_LC_heatmap

Rainbow_Emerald %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  Rainbow_Emerald_heatmap

Star_Emerald %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  Star_Emerald_heatmap

MacN_Emerald %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  MacN_Emerald_heatmap

Star_Rainbow %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  Star_Rainbow_heatmap

MacN_Rainbow %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  MacN_Rainbow_heatmap

MacN_Star %>%
  unite("gene_name", annot_ssid:annot_ofav, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  MacN_Star_heatmap


#### HEATMAPS TREATMENT ####

# only plotting comparisons with >15 shared, annotated DEGs
# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/ofav/host/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% CH_CC_heatmap$Protein_ofav)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% CH_CC_heatmap$Protein_ofav, CH_CC_heatmap$Orthogroup, CH_CC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/host/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% CH_CC_heatmap$Protein_ssid)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% CH_CC_heatmap$Protein_ssid, CH_CC_heatmap$Orthogroup, CH_CC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(CH_CC_heatmap$Orthogroup, CH_CC_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_CH_CC_p0.1.pdf", height=3.5, width=34)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(CH_CC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/ofav/host/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% LH_CC_heatmap$Protein_ofav)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% LH_CC_heatmap$Protein_ofav, LH_CC_heatmap$Orthogroup, LH_CC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/host/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% LH_CC_heatmap$Protein_ssid)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% LH_CC_heatmap$Protein_ssid, LH_CC_heatmap$Orthogroup, LH_CC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(LH_CC_heatmap$Orthogroup, LH_CC_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_LH_CC_p0.1.pdf", height=3.5, width=25)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(LH_CC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/ofav/host/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% CH_LC_heatmap$Protein_ofav)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% CH_LC_heatmap$Protein_ofav, CH_LC_heatmap$Orthogroup, CH_LC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/host/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% CH_LC_heatmap$Protein_ssid)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% CH_LC_heatmap$Protein_ssid, CH_LC_heatmap$Orthogroup, CH_LC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(CH_LC_heatmap$Orthogroup, CH_LC_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_CH_LC_p0.1.pdf", height=4.5, width=34)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(CH_LC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/ofav/host/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% LH_LC_heatmap$Protein_ofav)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% LH_LC_heatmap$Protein_ofav, LH_LC_heatmap$Orthogroup, LH_LC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/host/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% LH_LC_heatmap$Protein_ssid)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% LH_LC_heatmap$Protein_ssid, LH_LC_heatmap$Orthogroup, LH_LC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(LH_LC_heatmap$Orthogroup, LH_LC_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_LH_LC_p0.1.pdf", height=4.5, width=30)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(LH_LC_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()


#### HEATMAPS SITE ####

# only plotting comparisons with >15 shared, annotated DEGs
# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/ofav/host/vsd.RData")
design_ofav <- design
vsd_ofav <- subset(vsd, rownames(vsd) %in% MacN_Emerald_heatmap$Protein_ofav)

vsd_ofav %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ofav") %>%
  mutate(Protein_ofav = if_else(Protein_ofav %in% MacN_Emerald_heatmap$Protein_ofav, MacN_Emerald_heatmap$Orthogroup, MacN_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ofav") %>%
  as.matrix() -> vsd_ofav

load("../../DESeq2/ssid/host/vsd.RData")
design_ssid <- design
vsd_ssid <- subset(vsd, rownames(vsd) %in% MacN_Emerald_heatmap$Protein_ssid)

vsd_ssid %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_ssid") %>%
  mutate(Protein_ssid = if_else(Protein_ssid %in% MacN_Emerald_heatmap$Protein_ssid, MacN_Emerald_heatmap$Orthogroup, MacN_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_ssid") %>%
  as.matrix() -> vsd_ssid

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_ofav,vsd_ssid)
design_comb <- rbind(design_ofav,design_ssid)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(MacN_Emerald_heatmap$Orthogroup, MacN_Emerald_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_MacN_Emerald_p0.1.pdf", height=3.5, width=30)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(MacN_Emerald_heatmap$lpv_ssid)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_ofav)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# saving dataframes
save(orthologs, LC_CC, CH_CC, LH_CC, CH_LC, LH_CH, LH_LC, Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star, LC_CC_heatmap, CH_CC_heatmap, LH_CC_heatmap, CH_LC_heatmap, LH_CH_heatmap, LH_LC_heatmap, Rainbow_Emerald_heatmap, Star_Emerald_heatmap, MacN_Emerald_heatmap, Star_Rainbow_heatmap, MacN_Rainbow_heatmap, MacN_Star_heatmap, file = "orthofinder_DEGs.RData")


#### CHERRY PICKING ####

ofav_MacN_Emerald_cherry <- read.csv(file = "../../DESeq2/ofav/host/cherrypicking_MacN_Emerald.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_ofav_MacNEmerald"=lpv, "annot_ofav_MacNEmerald"=annot)

ofav_LH_CC_cherry <- read.csv(file = "../../DESeq2/ofav/host/cherrypicking_LH_CC.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_ofav_LHCC"=lpv, "annot_ofav_LHCC"=annot)

ssid_MacN_Emerald_cherry <- read.csv(file = "../../DESeq2/ssid/host/cherrypicking_MacN_Emerald.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_ssid_MacNEmerald"=lpv, "annot_ssid_MacNEmerald"=annot)

ssid_LH_CC_cherry <- read.csv(file = "../../DESeq2/ssid/host/cherrypicking_LH_CC.csv") %>%
  select(gene, lpv, annot) %>%
  rename("lpv_ssid_LHCC"=lpv, "annot_ssid_LHCC"=annot)

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  right_join(ofav_MacN_Emerald_cherry, by = c("Protein_ofav" = "gene")) %>%
  write.csv(file="ofav_MacN_Emerald_cherry.csv")

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  right_join(ofav_LH_CC_cherry, by = c("Protein_ofav" = "gene")) %>%
  write.csv(file="ofav_LH_CC_cherry.csv")

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  right_join(ssid_MacN_Emerald_cherry, by = c("Protein_ssid" = "gene")) %>%
  write.csv(file="ssid_MacN_Emerald_cherry.csv")

orthologs %>%
  rename("Protein_ofav" = 
           Ofaveolata_out_PRO, "Protein_ssid" = 	
           Siderastrea_out_PRO) %>%
  separate_rows(., Protein_ofav, sep = ",") %>%
  separate_rows(., Protein_ssid, sep = ",") %>%
  unique() %>%
  right_join(ssid_LH_CC_cherry, by = c("Protein_ssid" = "gene")) %>%
  write.csv(file="ssid_LH_CC_cherry.csv")


#### BOXPLOTS MATCHING ####

library(tidyverse)
library(stringr)
library(rcompanion)
library(rstatix)
library(ggpubr)
library(scales)

# DEGs with matching orthogroups
# OG0007874 (Spondin 2b, extracellular matrix protein)
Ofaveolata007753_site <- read.csv(file = "../../DESeq2/ofav/host/Ofaveolata007753_site.csv")
Ofaveolata007753_site$site <- factor(Ofaveolata007753_site$site, levels = c("Emerald", "Rainbow", "Star", "MacN"))

Ofaveolata007753_treat <- read.csv(file = "../../DESeq2/ofav/host/Ofaveolata007753_treat.csv")
Ofaveolata007753_treat$treat <- factor(Ofaveolata007753_treat$treat, levels = c("CC", "LC", "CH", "LH"))

Siderastrea444176_site <- read.csv(file = "../../DESeq2/ssid/host/Siderastrea444176_site.csv")
Siderastrea444176_site$site <- factor(Siderastrea444176_site$site, levels = c("Emerald", "Rainbow", "Star", "MacN"))

# ANOVA and Tukey's
Ofaveolata007753_site_stats <- aov(count~site,data=Ofaveolata007753_site) %>%
  tukey_hsd()
Ofaveolata007753_site_stats

Ofaveolata007753_treat_stats <- aov(count~treat,data=Ofaveolata007753_treat) %>%
  tukey_hsd()
Ofaveolata007753_treat_stats

Siderastrea444176_site_stats <- aov(count~site,data=Siderastrea444176_site) %>%
  tukey_hsd()
Siderastrea444176_site_stats

Ofaveolata007753_site_plot <-
  ggboxplot(
    Ofaveolata007753_site,
    x = "site",
    y = "count",
    color = "grey30",
    fill = "site",
    palette = c("#018571","#80cdc1","#dfc27d","#a6611a"),
    title = "O. faveolata site",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'site') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") +
  scale_y_log10(limit = c(0.25,100000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata007753_site_stats,label="p.adj.signif",y.position=3.75,step.increase=0.075,inherit.aes=FALSE,size=3)
Ofaveolata007753_site_plot

Ofaveolata007753_treat_plot <-
  ggboxplot(
    Ofaveolata007753_treat,
    x = "treat",
    y = "count",
    color = "grey30",
    fill = "treat",
    palette = c("#92c5de","#f4a582","#0571b0","#ca0020"),
    title = "O. faveolata treat",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'treat') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") +
  scale_y_log10(limit = c(0.25,100000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Ofaveolata007753_treat_stats,label="p.adj.signif",y.position=3.75,step.increase=0.075,inherit.aes=FALSE,size=3)
Ofaveolata007753_treat_plot

Siderastrea444176_site_plot <-
  ggboxplot(
    Siderastrea444176_site,
    x = "site",
    y = "count",
    color = "grey30",
    fill = "site",
    palette = c("#018571","#80cdc1","#dfc27d","#a6611a"),
    title = "S. siderea site",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = element_blank(),
           y = "Log Normalized Counts",
           fill = 'site') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "none") +
  scale_y_log10(limit = c(0.25,100000), breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  stat_pvalue_manual(Siderastrea444176_site_stats,label="p.adj.signif",y.position=3,step.increase=0.1,inherit.aes=FALSE,size=3)
Siderastrea444176_site_plot


#### MULTIPLOT MATCHING ####

species_match<-ggarrange(Ofaveolata007753_site_plot,
                         Ofaveolata007753_treat_plot,
                         Siderastrea444176_site_plot,
                              heights = c(2),
                              widths = c(6,6,6),
                              ncol = 3,
                              nrow = 1)
species_match<-annotate_figure(species_match, top = text_grob("Spondin 2b, extracellular matrix protein", color = "black", face = "bold", size = 14))
species_match

ggsave("species_match.pdf", species_match, width=12, height=4,dpi = 300)
