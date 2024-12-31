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
# Check the number of genes assigned to orthogroups (e.g., 'OrthoFinder assigned 25494 genes (90.2% of total) to 7908 orthogroups')
# Ideally, it should be >80%


#### ORTHOLOGS ####

# load("orthofinder_DEGs.RData") # if previously run
orthologs <- read.table(file = "orthofinder/OrthoFinder/Results_Oct09/Orthologues/Orthologues_Durusdinium_out_PRO/Durusdinium_out_PRO__v__Breviolum_out_PRO.tsv", sep = "\t", header = TRUE, quote="", fill=FALSE)

orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  mutate(Protein_symD = trimws(Protein_symD)) %>%
  mutate(Protein_symB = trimws(Protein_symB)) %>%
  unique() -> orthologs_unique
write.table(orthologs_unique, file = "orthologs_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### IMPORT RAW COUNTS ####

ofav_symD <- read.table(file = "../../../raw/ofav/allcounts_sym.txt") %>%
  tibble::rownames_to_column("gene")

ofav_symB <- read.table(file = "../../../raw/ofav/allcounts_sym2.txt") %>%
  tibble::rownames_to_column("gene")


#### DEG MATCHING TREATMENT ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds symD and symB gene annotations, and 5) then pulls on corresponding KOG classes

# LC vs CC for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(ofav_symD, by = c("Protein_symD" = "gene")) %>%
  inner_join(ofav_symB, by = c("Protein_symB" = "gene")) -> orthologs_ofav

%>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="LC_CC", .before="Orthogroup") -> LC_CC
LC_CC$Orthogroup <- make.unique(LC_CC$Orthogroup, sep = "_") 


#### DESEQ IMPORT TREATMENT ####

symD_LC_CC_lpv <- read.csv(file = "../../DESeq2/symD/sym/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_CH_CC_lpv <- read.csv(file = "../../DESeq2/symD/sym/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_LH_CC_lpv <- read.csv(file = "../../DESeq2/symD/sym/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_CH_LC_lpv <- read.csv(file = "../../DESeq2/symD/sym/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_LH_CH_lpv <- read.csv(file = "../../DESeq2/symD/sym/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_LH_LC_lpv <- read.csv(file = "../../DESeq2/symD/sym/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symB_LC_CC_lpv <- read.csv(file = "../../DESeq2/symD/sym2/LC_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_CH_CC_lpv <- read.csv(file = "../../DESeq2/symD/sym2/CH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_LH_CC_lpv <- read.csv(file = "../../DESeq2/symD/sym2/LH_CC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_CH_LC_lpv <- read.csv(file = "../../DESeq2/symD/sym2/CH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_LH_CH_lpv <- read.csv(file = "../../DESeq2/symD/sym2/LH_CH_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_LH_LC_lpv <- read.csv(file = "../../DESeq2/symD/sym2/LH_LC_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)


#### DESEQ IMPORT SITE ####
symD_Rainbow_Emerald_lpv <- read.csv(file = "../../DESeq2/symD/sym/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_Star_Emerald_lpv <- read.csv(file = "../../DESeq2/symD/sym/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_MacN_Emerald_lpv <- read.csv(file = "../../DESeq2/symD/sym/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_Star_Rainbow_lpv <- read.csv(file = "../../DESeq2/symD/sym/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_MacN_Rainbow_lpv <- read.csv(file = "../../DESeq2/symD/sym/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symD_MacN_Star_lpv <- read.csv(file = "../../DESeq2/symD/sym/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symD" = lpv)

symB_Rainbow_Emerald_lpv <- read.csv(file = "../../DESeq2/symD/sym2/Rainbow_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_Star_Emerald_lpv <- read.csv(file = "../../DESeq2/symD/sym2/Star_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_MacN_Emerald_lpv <- read.csv(file = "../../DESeq2/symD/sym2/MacN_Emerald_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_Star_Rainbow_lpv <- read.csv(file = "../../DESeq2/symD/sym2/Star_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_MacN_Rainbow_lpv <- read.csv(file = "../../DESeq2/symD/sym2/MacN_Rainbow_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)

symB_MacN_Star_lpv <- read.csv(file = "../../DESeq2/symD/sym2/MacN_Star_lpv.csv") %>%
  select(gene, lpv) %>%
  rename("lpv_symB" = lpv)


#### DEG MATCHING TREATMENT ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds symD and symB gene annotations, and 5) then pulls on corresponding KOG classes

# LC vs CC for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_LC_CC_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_LC_CC_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="LC_CC", .before="Orthogroup") -> LC_CC
LC_CC$Orthogroup <- make.unique(LC_CC$Orthogroup, sep = "_") 

# CH vs CC for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_CH_CC_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_CH_CC_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="CH_CC", .before="Orthogroup") -> CH_CC
CH_CC$Orthogroup <- make.unique(CH_CC$Orthogroup, sep = "_") 

# LH vs CC for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_LH_CC_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_LH_CC_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="LH_CC", .before="Orthogroup") -> LH_CC
LH_CC$Orthogroup <- make.unique(LH_CC$Orthogroup, sep = "_") 

# CH vs LC for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_CH_LC_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_CH_LC_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="CH_LC", .before="Orthogroup") -> CH_LC
CH_LC$Orthogroup <- make.unique(CH_LC$Orthogroup, sep = "_") 

# LH vs CH for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_LH_CH_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_LH_CH_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="LH_CH", .before="Orthogroup") -> LH_CH
LH_CH$Orthogroup <- make.unique(LH_CH$Orthogroup, sep = "_") 

# LH vs LC for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_LH_LC_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_LH_LC_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="LH_LC", .before="Orthogroup") -> LH_LC
LH_LC$Orthogroup <- make.unique(LH_LC$Orthogroup, sep = "_") 

# joining all matching DEGs into a single dataframe
orthofinder_treatment <- bind_rows(LC_CC,CH_CC,LH_CC,CH_LC,LH_CH,LH_LC)
write.csv(orthofinder_treatment, file="orthofinder_treatment.csv")


#### DEG MATCHING SITE ####

# This section of code does several things: 1) rename common orthologs, 2) join with -log10(pval), 3) filter by 0.1 pval cutoff (log10(0.1)=1), 4) adds symD and symB gene annotations, and 5) then pulls on corresponding KOG classes

# Rainbow vs Emerald for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_Rainbow_Emerald_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_Rainbow_Emerald_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="Rainbow_Emerald", .before="Orthogroup") -> Rainbow_Emerald
Rainbow_Emerald$Orthogroup <- make.unique(Rainbow_Emerald$Orthogroup, sep = "_") 

# Star vs Emerald for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_Star_Emerald_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_Star_Emerald_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="Star_Emerald", .before="Orthogroup") -> Star_Emerald
Star_Emerald$Orthogroup <- make.unique(Star_Emerald$Orthogroup, sep = "_") 

# MacN vs Emerald for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_MacN_Emerald_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_MacN_Emerald_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="MacN_Emerald", .before="Orthogroup") -> MacN_Emerald
MacN_Emerald$Orthogroup <- make.unique(MacN_Emerald$Orthogroup, sep = "_") 

# Star vs Rainbow for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_Star_Rainbow_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_Star_Rainbow_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="Star_Rainbow", .before="Orthogroup") -> Star_Rainbow
Star_Rainbow$Orthogroup <- make.unique(Star_Rainbow$Orthogroup, sep = "_") 

# MacN vs Rainbow for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_MacN_Rainbow_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_MacN_Rainbow_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="MacN_Rainbow", .before="Orthogroup") -> MacN_Rainbow
MacN_Rainbow$Orthogroup <- make.unique(MacN_Rainbow$Orthogroup, sep = "_") 

# MacN vs Star for both species
orthologs %>%
  rename("Protein_symD" = 
           Durusdinium_out_PRO, "Protein_symB" = 	
           Breviolum_out_PRO) %>%
  separate_rows(., Protein_symD, sep = ",") %>%
  separate_rows(., Protein_symB, sep = ",") %>%
  unique() %>%
  inner_join(symB_MacN_Star_lpv, by = c("Protein_symB" = "gene")) %>%
  inner_join(symD_MacN_Star_lpv, by = c("Protein_symD" = "gene")) %>%
  filter(abs(lpv_symB) >= 1 & abs(lpv_symD) >= 1) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symD/shoguchi/Durusdinium_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symD = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symD" = "gene")) %>%
  left_join(read.table(file = "../../../Annotations/symB/magana2/Breviolum_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG_symB = V2) %>%
              dplyr::select(-V1, -V2), by = c("Protein_symB" = "gene")) %>%
  mutate(comparison="MacN_Star", .before="Orthogroup") -> MacN_Star
MacN_Star$Orthogroup <- make.unique(MacN_Star$Orthogroup, sep = "_") 

# joining all matching DEGs into a single dataframe
orthofinder_site <- bind_rows(Rainbow_Emerald,Star_Emerald,MacN_Emerald,Star_Rainbow,MacN_Rainbow,MacN_Star)
write.csv(orthofinder_site, file="orthofinder_site.csv")


#### VENN DIAGRAMS TREATMENT ####

# first creating a set of up/downregulated DEGs by species
CH_CC %>%
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

CH_CC %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

CH_CC %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

CH_CC %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

LH_CC %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

LH_CC %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

LH_CC %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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


#### VENN DIAGRAMS SITE ####

# first creating a set of up/downregulated DEGs by species
Rainbow_Emerald %>%
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

Rainbow_Emerald %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

Rainbow_Emerald %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

Rainbow_Emerald %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

Star_Emerald %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

Star_Emerald %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

Star_Emerald %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

MacN_Emerald %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

MacN_Emerald %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

MacN_Emerald %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

Star_Rainbow %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

Star_Rainbow %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

Star_Rainbow %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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
  filter(lpv_symB >= 1) %>%
  pull(Orthogroup) -> symB_up

MacN_Rainbow %>%
  filter(lpv_symB <= -1) %>%
  pull(Orthogroup) -> symB_down

MacN_Rainbow %>%
  filter(lpv_symD >= 1) %>%
  pull(Orthogroup) -> symD_up

MacN_Rainbow %>%
  filter(lpv_symD <= -1) %>%
  pull(Orthogroup) -> symD_down

venn=venn.diagram(
  x = list("B up"=symB_up, "B down"=symB_down,"D up"=symD_up, "D down"=symD_down),
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
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  CH_CC_heatmap

LH_CC %>%
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  LH_CC_heatmap

Rainbow_Emerald %>%
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  Rainbow_Emerald_heatmap

Star_Emerald %>%
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  Star_Emerald_heatmap

MacN_Emerald %>%
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  MacN_Emerald_heatmap

Star_Rainbow %>%
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  Star_Rainbow_heatmap

MacN_Rainbow %>%
  unite("gene_name", annot_symD:annot_symB, sep = " / ", remove = FALSE) %>%
  mutate(gene_name = str_replace(gene_name, "NA / NA","")) %>%
  mutate(gene_name = str_replace(gene_name, "- / -","")) %>%
  mutate(gene_name = na_if(gene_name,"")) %>%
  filter(!is.na(gene_name)) ->  MacN_Rainbow_heatmap


#### HEATMAPS TREATMENT ####

# only plotting comparisons with >10 shared, annotated DEGs
# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/symD/sym/vsd.RData")
design_symD <- design
vsd_symD <- subset(vsd, rownames(vsd) %in% LH_CC_heatmap$Protein_symD)

vsd_symD %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symD") %>%
  mutate(Protein_symD = if_else(Protein_symD %in% LH_CC_heatmap$Protein_symD, LH_CC_heatmap$Orthogroup, LH_CC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symD") %>%
  as.matrix() -> vsd_symD

load("../../DESeq2/symD/sym2/vsd.RData")
design_symB <- design
vsd_symB <- subset(vsd, rownames(vsd) %in% LH_CC_heatmap$Protein_symB)

vsd_symB %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symB") %>%
  mutate(Protein_symB = if_else(Protein_symB %in% LH_CC_heatmap$Protein_symB, LH_CC_heatmap$Orthogroup, LH_CC_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symB") %>%
  as.matrix() -> vsd_symB

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_symD,vsd_symB)
design_comb <- rbind(design_symD,design_symB)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(LH_CC_heatmap$Orthogroup, LH_CC_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_LH_CC_p0.1.pdf", height=2.5, width=30)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(LH_CC_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()


#### HEATMAPS SITE ####

# only plotting comparisons with >10 shared, annotated DEGs
# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/symD/sym/vsd.RData")
design_symD <- design
vsd_symD <- subset(vsd, rownames(vsd) %in% Rainbow_Emerald_heatmap$Protein_symD)

vsd_symD %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symD") %>%
  mutate(Protein_symD = if_else(Protein_symD %in% Rainbow_Emerald_heatmap$Protein_symD, Rainbow_Emerald_heatmap$Orthogroup, Rainbow_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symD") %>%
  as.matrix() -> vsd_symD

load("../../DESeq2/symD/sym2/vsd.RData")
design_symB <- design
vsd_symB <- subset(vsd, rownames(vsd) %in% Rainbow_Emerald_heatmap$Protein_symB)

vsd_symB %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symB") %>%
  mutate(Protein_symB = if_else(Protein_symB %in% Rainbow_Emerald_heatmap$Protein_symB, Rainbow_Emerald_heatmap$Orthogroup, Rainbow_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symB") %>%
  as.matrix() -> vsd_symB

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_symD,vsd_symB)
design_comb <- rbind(design_symD,design_symB)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(Rainbow_Emerald_heatmap$Orthogroup, Rainbow_Emerald_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.01 
pdf(file="heatmap_Rainbow_Emerald_p0.01.pdf", height=45, width=75)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(Rainbow_Emerald_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-2, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# p < 0.001
pdf(file="heatmap_Rainbow_Emerald_p0.001.pdf", height=2.5, width=45)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(Rainbow_Emerald_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-3, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off()

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/symD/sym/vsd.RData")
design_symD <- design
vsd_symD <- subset(vsd, rownames(vsd) %in% Star_Emerald_heatmap$Protein_symD)

vsd_symD %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symD") %>%
  mutate(Protein_symD = if_else(Protein_symD %in% Star_Emerald_heatmap$Protein_symD, Star_Emerald_heatmap$Orthogroup, Star_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symD") %>%
  as.matrix() -> vsd_symD

load("../../DESeq2/symD/sym2/vsd.RData")
design_symB <- design
vsd_symB <- subset(vsd, rownames(vsd) %in% Star_Emerald_heatmap$Protein_symB)

vsd_symB %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symB") %>%
  mutate(Protein_symB = if_else(Protein_symB %in% Star_Emerald_heatmap$Protein_symB, Star_Emerald_heatmap$Orthogroup, Star_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symB") %>%
  as.matrix() -> vsd_symB

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_symD,vsd_symB)
design_comb <- rbind(design_symD,design_symB)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(Star_Emerald_heatmap$Orthogroup, Star_Emerald_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.0001
pdf(file="heatmap_Star_Emerald_p1e-12.pdf", height=2.5, width=45)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(Star_Emerald_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-12, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off() # too many genes

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/symD/sym/vsd.RData")
design_symD <- design
vsd_symD <- subset(vsd, rownames(vsd) %in% MacN_Emerald_heatmap$Protein_symD)

vsd_symD %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symD") %>%
  mutate(Protein_symD = if_else(Protein_symD %in% MacN_Emerald_heatmap$Protein_symD, MacN_Emerald_heatmap$Orthogroup, MacN_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symD") %>%
  as.matrix() -> vsd_symD

load("../../DESeq2/symD/sym2/vsd.RData")
design_symB <- design
vsd_symB <- subset(vsd, rownames(vsd) %in% MacN_Emerald_heatmap$Protein_symB)

vsd_symB %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symB") %>%
  mutate(Protein_symB = if_else(Protein_symB %in% MacN_Emerald_heatmap$Protein_symB, MacN_Emerald_heatmap$Orthogroup, MacN_Emerald_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symB") %>%
  as.matrix() -> vsd_symB

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_symD,vsd_symB)
design_comb <- rbind(design_symD,design_symB)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(MacN_Emerald_heatmap$Orthogroup, MacN_Emerald_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_MacN_Emerald_p1e-6.pdf", height=3.5, width=30)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(MacN_Emerald_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-6, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off() # too many genes

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/symD/sym/vsd.RData")
design_symD <- design
vsd_symD <- subset(vsd, rownames(vsd) %in% Star_Rainbow_heatmap$Protein_symD)

vsd_symD %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symD") %>%
  mutate(Protein_symD = if_else(Protein_symD %in% Star_Rainbow_heatmap$Protein_symD, Star_Rainbow_heatmap$Orthogroup, Star_Rainbow_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symD") %>%
  as.matrix() -> vsd_symD

load("../../DESeq2/symD/sym2/vsd.RData")
design_symB <- design
vsd_symB <- subset(vsd, rownames(vsd) %in% Star_Rainbow_heatmap$Protein_symB)

vsd_symB %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symB") %>%
  mutate(Protein_symB = if_else(Protein_symB %in% Star_Rainbow_heatmap$Protein_symB, Star_Rainbow_heatmap$Orthogroup, Star_Rainbow_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symB") %>%
  as.matrix() -> vsd_symB

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_symD,vsd_symB)
design_comb <- rbind(design_symD,design_symB)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(Star_Rainbow_heatmap$Orthogroup, Star_Rainbow_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_Star_Rainbow_p0.1.pdf", height=8, width=48)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(Star_Rainbow_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off() 

# first loading variance stabilized arrays of gene counts, then replacing species-specific gene IDs with orthogroup ID
load("../../DESeq2/symD/sym/vsd.RData")
design_symD <- design
vsd_symD <- subset(vsd, rownames(vsd) %in% MacN_Rainbow_heatmap$Protein_symD)

vsd_symD %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symD") %>%
  mutate(Protein_symD = if_else(Protein_symD %in% MacN_Rainbow_heatmap$Protein_symD, MacN_Rainbow_heatmap$Orthogroup, MacN_Rainbow_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symD") %>%
  as.matrix() -> vsd_symD

load("../../DESeq2/symD/sym2/vsd.RData")
design_symB <- design
vsd_symB <- subset(vsd, rownames(vsd) %in% MacN_Rainbow_heatmap$Protein_symB)

vsd_symB %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein_symB") %>%
  mutate(Protein_symB = if_else(Protein_symB %in% MacN_Rainbow_heatmap$Protein_symB, MacN_Rainbow_heatmap$Orthogroup, MacN_Rainbow_heatmap$Orthogroup)) %>%
  column_to_rownames(var = "Protein_symB") %>%
  as.matrix() -> vsd_symB

# combining both matrices and design metadata for plotting
vsd_comb <- cbind(vsd_symD,vsd_symB)
design_comb <- rbind(design_symD,design_symB)
design_comb$id <- as.factor(gsub("-",".", design_comb$id))
design_comb$full_id <- paste(design_comb$id,design_comb$site,design_comb$treat,sep=".")

# Make sure the 'uniHeatmap.R' script is in your working directory
source("uniHeatmap.R")

# creating a lookup table of orthogroup to gene annotations
gene_names <- as.data.frame(cbind(MacN_Rainbow_heatmap$Orthogroup, MacN_Rainbow_heatmap$gene_name))

# heatmaps
# cutoff -1 (0.1), -1.3 (0.05), -2 (0.01), -3 (0.001), -6 (1e6)
# p < 0.1 (all genes)
pdf(file="heatmap_MacN_Rainbow_p0.1.pdf", height=7, width=48)
uniHeatmap(vsd=vsd_comb,gene.names=gene_names,
           metric=-(abs(MacN_Rainbow_heatmap$lpv_symB)), # metric of gene significance
           # metric2=-(abs(MacN_Emerald$lpv_symD)),
           cutoff=-1, 
           sort=c(1:ncol(vsd_comb)), # overrides sorting of columns according to hierarchical clustering
           # sort=order(design_comb$full_id), 
           cex=0.8,
           pdf=F,
)
dev.off() 

# saving dataframes
save(orthologs, LC_CC, CH_CC, LH_CC, CH_LC, LH_CH, LH_LC, Rainbow_Emerald, Star_Emerald, MacN_Emerald, Star_Rainbow, MacN_Rainbow, MacN_Star, LC_CC_heatmap, CH_CC_heatmap, LH_CC_heatmap, Rainbow_Emerald_heatmap, Star_Emerald_heatmap, MacN_Emerald_heatmap, Star_Rainbow_heatmap, MacN_Rainbow_heatmap, file = "orthofinder_DEGs.RData")


#### KOG MATCHING SITE ####

# Rainbow vs Emerald for both species
Rainbow_Emerald %>%
  mutate(KOG_symD = replace(KOG_symD, KOG_symD == "", NA)) %>%
  filter(lpv_symD >= 1) %>%
  count(KOG_symD) %>%
  rename("KOG" = KOG_symD, "symD_up" = n) -> KOG_symD_up

Rainbow_Emerald %>%
  mutate(KOG_symD = replace(KOG_symD, KOG_symD == "", NA)) %>%
  filter(lpv_symD <= -1) %>%
  count(KOG_symD) %>%
  rename("KOG" = KOG_symD, "symD_down" = n) -> KOG_symD_down

Rainbow_Emerald %>%
  mutate(KOG_symB = replace(KOG_symB, KOG_symB == "", NA)) %>%
  filter(lpv_symB >= 1) %>%
  count(KOG_symB) %>%
  rename("KOG" = KOG_symB, "symB_up" = n) -> KOG_symB_up

Rainbow_Emerald %>%
  mutate(KOG_symB = replace(KOG_symB, KOG_symB == "", NA)) %>%
  filter(lpv_symB <= -1) %>%
  count(KOG_symB) %>%
  rename("KOG" = KOG_symB, "symB_down" = n) -> KOG_symB_down

# joining all KOG class sums in a single dataframe
KOG_symD_up %>%
  left_join(KOG_symB_up, by = "KOG") %>%
  inner_join(KOG_symD_down, by = "KOG") %>%
  inner_join(KOG_symB_down, by = "KOG") -> KOG_match

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
  mutate(KOG_symD = replace(KOG_symD, KOG_symD == "", NA)) %>%
  filter(lpv_symD >= 1) %>%
  count(KOG_symD) %>%
  rename("KOG" = KOG_symD, "symD_up" = n) -> KOG_symD_up

Star_Emerald %>%
  mutate(KOG_symD = replace(KOG_symD, KOG_symD == "", NA)) %>%
  filter(lpv_symD <= -1) %>%
  count(KOG_symD) %>%
  rename("KOG" = KOG_symD, "symD_down" = n) -> KOG_symD_down

Star_Emerald %>%
  mutate(KOG_symB = replace(KOG_symB, KOG_symB == "", NA)) %>%
  filter(lpv_symB >= 1) %>%
  count(KOG_symB) %>%
  rename("KOG" = KOG_symB, "symB_up" = n) -> KOG_symB_up

Star_Emerald %>%
  mutate(KOG_symB = replace(KOG_symB, KOG_symB == "", NA)) %>%
  filter(lpv_symB <= -1) %>%
  count(KOG_symB) %>%
  rename("KOG" = KOG_symB, "symB_down" = n) -> KOG_symB_down

# joining all KOG class sums in a single dataframe
KOG_symD_up %>%
  left_join(KOG_symB_up, by = "KOG") %>%
  inner_join(KOG_symD_down, by = "KOG") %>%
  inner_join(KOG_symB_down, by = "KOG") -> KOG_match

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
  mutate(KOG_symD = replace(KOG_symD, KOG_symD == "", NA)) %>%
  filter(lpv_symD >= 1) %>%
  count(KOG_symD) %>%
  rename("KOG" = KOG_symD, "symD_up" = n) -> KOG_symD_up

MacN_Emerald %>%
  mutate(KOG_symD = replace(KOG_symD, KOG_symD == "", NA)) %>%
  filter(lpv_symD <= -1) %>%
  count(KOG_symD) %>%
  rename("KOG" = KOG_symD, "symD_down" = n) -> KOG_symD_down

MacN_Emerald %>%
  mutate(KOG_symB = replace(KOG_symB, KOG_symB == "", NA)) %>%
  filter(lpv_symB >= 1) %>%
  count(KOG_symB) %>%
  rename("KOG" = KOG_symB, "symB_up" = n) -> KOG_symB_up

MacN_Emerald %>%
  mutate(KOG_symB = replace(KOG_symB, KOG_symB == "", NA)) %>%
  filter(lpv_symB <= -1) %>%
  count(KOG_symB) %>%
  rename("KOG" = KOG_symB, "symB_down" = n) -> KOG_symB_down

# joining all KOG class sums in a single dataframe
KOG_symD_up %>%
  left_join(KOG_symB_up, by = "KOG") %>%
  inner_join(KOG_symD_down, by = "KOG") %>%
  inner_join(KOG_symB_down, by = "KOG") -> KOG_match

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

# too few matching genes for MacN vs Rainbow in both species

# no matching genes for MacN vs Star in both species
