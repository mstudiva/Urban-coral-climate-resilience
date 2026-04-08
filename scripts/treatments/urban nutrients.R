#### packages ####

library(tidyverse)
library(dplyr)
library(ggpubr)
library(vegan)
library(ape)


#### temp data import ####

nutrients <- read.csv(file="../../data/treatments/urban nutrients.csv", head=T) 
nutrients$Tank<- as.factor(nutrients$Tank) 
nutrients$Time <- as.factor(nutrients$Time)
nutrients$Species <- as.factor(nutrients$Species)
nutrients$Site <- as.factor(nutrients$Site)
nutrients$Treatment <- as.factor(nutrients$Treatment)
str(nutrients)

nutrients %>%
  group_by(Species,Time, Treatment) %>%
  summarize(meanSi=mean(Si), 
            lowSi = min(Si), 
            highSi= max(Si),
            meanNO2=mean(NO2), 
            lowNO2 = min(NO2), 
            highNO2= max(NO2),
            meanNO3=mean(NO3), 
            lowNO3 = min(NO3), 
            highNO3= max(NO3),
            meanNH4=mean(NH4), 
            lowNH4 = min(NH4), 
            highNH4= max(NH4),
            meanPO4=mean(PO4), 
            lowPO4 = min(PO4), 
            highPO4= max(PO4))
            