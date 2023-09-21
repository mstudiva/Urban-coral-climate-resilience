#### packages ####

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(rstatix)
library(MASS)


#### data import ####

ipam <- read.csv(file="urban ipam.csv", head=T)
str(ipam)

# changing genotype numbers to be distinct across sites/species
genotypes <- read.csv(file="genotype lookup table.csv", head=T)
genotypes

ipam %>%
  inner_join(genotypes, by = "sampleID") -> ipam

# coding factors
ipam$geno <- as.factor(ipam$geno)
ipam$tank <- as.factor(ipam$tank)
ipam$time=factor(ipam$time, levels=c("Initial", "0", "3", "7", "10", "14", "17", "21", "24", "28", "31", "35", "39", "42", "45")) 
ipam$site=factor(ipam$site, levels=c("Emerald", "Rainbow", "Star", "MacN")) 
ipam$treatment=factor(ipam$treatment, levels=c("CC", "CH", "LC", "LH")) 
str(ipam)


#### data normality and assumption testing ####

# subsetting data frames by species
ipam_ofav <- subset(ipam, species=="Ofav")
ipam_ssid <- subset(ipam, species=="Ssid")
str(ipam_ofav)
str(ipam_ssid)

# now removing the initial data points for statistical tests, since they are all '1'
ipam_ofav_sub <- subset(ipam_ofav, ratio!="1")
ipam_ssid_sub <- subset(ipam_ssid, ratio!="1")

# recoding time as numeric for plotting/testing
ipam_ofav_sub$time  = as.numeric(ipam_ofav_sub$time)
ipam_ssid_sub$time  = as.numeric(ipam_ssid_sub$time)
str(ipam_ofav_sub)
str(ipam_ssid_sub)

# scatterplots of measurements through time for factors
pdf("urban ipam ofav normality.pdf")
ggscatter(
  ipam_ofav_sub, x = "time", y = "ratio",
  facet.by  = c("treatment", "site")) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.9)
dev.off()

pdf("urban ipam ssid normality.pdf")
ggscatter(
  ipam_ssid_sub, x = "time", y = "ratio",
  facet.by  = c("treatment", "site")) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.9)
dev.off()
# somewhat linear decline in all combinations of site + treatment over time

# testing for homogeneity of regression slopes
# check that interaction terms with covariable (time) are not significant
ipam_ofav_sub %>%
  anova_test(
    ratio ~ time + treatment + site + 
      treatment*site + time*treatment +
      time*site + time*site*treatment 
  )
# all time interactions nonsignificant

ipam_ssid_sub %>%
  anova_test(
    ratio ~ time + treatment + site + 
      treatment*site + time*treatment +
      time*site + time*site*treatment
  )
# time:treatment interaction significant

# normality of residuals
ofav_model <- lm(ratio ~ time + treatment*site, data = ipam_ofav_sub)
ofav_model_metrics <- augment(ofav_model)

ssid_model <- lm(ratio ~ time + treatment*site, data = ipam_ssid_sub)
ssid_model_metrics <- augment(ssid_model)

# assess normality of residuals using Shapiro-Wilk test
shapiro_test(ofav_model_metrics$.resid)
shapiro_test(ssid_model_metrics$.resid)
# both significant

# assess homogeneity of variances using Levene's Test
levene_test(.resid ~ treatment*site, data = ofav_model_metrics)
levene_test(.resid ~ treatment*site, data = ssid_model_metrics)
# both significant

# identifying outliers
ofav_model_metrics %>% 
  filter(abs(.std.resid) > 3) -> ofav_outliers
ssid_model_metrics %>% 
  filter(abs(.std.resid) > 3) -> ssid_outliers

# remove outliers if needed
ipam_ofav_sub <- ipam_ofav_sub[!rownames(ipam_ofav_sub) %in% ofav_outliers$.rownames,]
ipam_ssid_sub <- ipam_ssid_sub[!rownames(ipam_ssid_sub) %in% ssid_outliers$.rownames,]

#### ANCOVA ####

ofav_aov <- ipam_ofav_sub %>% 
  anova_test(ratio ~ time + treatment*site)
get_anova_table(ofav_aov) # all factors and treatment:site interaction significant
capture.output(get_anova_table(ofav_aov), file = "urban ipam ofav ancova.txt")

ssid_aov <- ipam_ssid_sub %>% 
  anova_test(ratio ~ time + treatment*site)
get_anova_table(ssid_aov) # all factors and treatment:site interaction significant
capture.output(get_anova_table(ssid_aov), file = "urban ipam ssid ancova.txt")


#### subsetting for controls ####

ipam_ofav %>%
  filter(time!="Initial") %>%
  mutate(time = as.character(time)) %>%
  mutate(time = as.numeric(time)) %>%
  filter(time<21) -> ipam_ofav_control

ipam_ssid %>%
  filter(time!="Initial") %>%
  mutate(time = as.character(time)) %>%
  mutate(time = as.numeric(time)) %>%
  filter(time<21) -> ipam_ssid_control

ggscatter(
  ipam_ofav_control, x = "time", y = "ratio",
  facet.by  = c("treatment", "site")) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.9)

ggscatter(
  ipam_ofav_control, x = "time", y = "ratio",
  facet.by  = c("treatment", "site")) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.9)
# somewhat linear decline in all combinations of site + treatment over time

# testing for homogeneity of regression slopes
# check that interaction terms with covariable (time) are not signficant
ipam_ofav_control %>%
  anova_test(
    ratio ~ time + treatment + site + 
      treatment*site + time*treatment +
      time*site + time*site*treatment
  )
# all time interactions nonsignificant

ipam_ssid_control %>%
  anova_test(
    ratio ~ time + treatment + site + 
      treatment*site + time*treatment +
      time*site + time*site*treatment
  )
# time:treatment interaction significant

# normality of residuals
ofav_model_control <- lm(ratio ~ time + treatment*site, data = ipam_ofav_control)
ofav_model_control_metrics <- augment(ofav_model_control)

ssid_model_control <- lm(ratio ~ time + treatment*site, data = ipam_ssid_control)
ssid_model_control_metrics <- augment(ssid_model_control) 

# assess normality of residuals using Shapiro-Wilk test
shapiro_test(ofav_model_control_metrics$.resid)
shapiro_test(ssid_model_control_metrics$.resid)
# both significant

# assess homogeneity of variances using Levene's Test
levene_test(.resid ~ treatment*site, data = ofav_model_control_metrics)
levene_test(.resid ~ treatment*site, data = ssid_model_control_metrics)
# both significant

# identifying outliers
ofav_model_control_metrics %>% 
  filter(abs(.std.resid) > 3) -> ofav_control_outliers
ssid_model_control_metrics %>% 
  filter(abs(.std.resid) > 3) -> ssid_control_outliers

# remove outliers
ipam_ofav_control <- ipam_ofav_control[!rownames(ipam_ofav_control) %in% ofav_control_outliers$.rownames,]
ipam_ssid_control <- ipam_ssid_control[!rownames(ipam_ssid_control) %in% ssid_control_outliers$.rownames,]

#### ANCOVA with controls ####

ofav_control_aov <- ipam_ofav_control %>% 
  anova_test(ratio ~ time + treatment*site)
get_anova_table(ofav_control_aov) # all factors and treatment:site interaction significant
capture.output(get_anova_table(ofav_control_aov), file = "urban ipam ofav controls ancova.txt")

ssid_control_aov <- ipam_ssid_control %>% 
  anova_test(ratio ~ time + treatment*site)
get_anova_table(ssid_control_aov) # time and treatment significant
capture.output(get_anova_table(ssid_control_aov), file = "urban ipam ssid controls ancova.txt")


#### subsetting for heat ####

ipam_ofav %>%
  filter(time!="Initial") %>%
  mutate(time = as.character(time)) %>%
  mutate(time = as.numeric(time)) %>%
  filter(temp!="control temp") -> ipam_ofav_heat

ipam_ssid %>%
  filter(time!="Initial") %>%
  mutate(time = as.character(time)) %>%
  mutate(time = as.numeric(time)) %>%
  filter(temp!="control temp") -> ipam_ssid_heat

ggscatter(
  ipam_ofav_heat, x = "time", y = "ratio",
  facet.by  = c("treatment", "site")) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.9)

ggscatter(
  ipam_ofav_heat, x = "time", y = "ratio",
  facet.by  = c("treatment", "site")) +
  geom_smooth(formula = y ~ x, method = "loess", span = 0.9)
# somewhat linear decline in all combinations of site + treatment over time

# testing for homogeneity of regression slopes
# check that interaction terms with covariable (time) are not signficant
ipam_ofav_heat %>%
  anova_test(
    ratio ~ time + treatment + site + 
      treatment*site + time*treatment +
      time*site + time*site*treatment
  )
# time:treatment interaction significant

ipam_ssid_heat %>%
  anova_test(
    ratio ~ time + treatment + site + 
      treatment*site + time*treatment +
      time*site + time*site*treatment
  )
# time:treatment interaction significant

# normality of residuals
ofav_model_heat <- lm(ratio ~ time + treatment*site, data = ipam_ofav_heat)
ofav_model_heat_metrics <- augment(ofav_model_heat)

ssid_model_heat <- lm(ratio ~ time + treatment*site, data = ipam_ssid_heat)
ssid_model_heat_metrics <- augment(ssid_model_heat)

# assess normality of residuals using Shapiro-Wilk test
shapiro_test(ofav_model_heat_metrics$.resid)
shapiro_test(ssid_model_heat_metrics$.resid)
# both significant

# assess homogeneity of variances using Levene's Test
levene_test(.resid ~ treatment*site, data = ofav_model_heat_metrics)
levene_test(.resid ~ treatment*site, data = ssid_model_heat_metrics)
# both significant

# identifying outliers
ofav_model_heat_metrics %>% 
  filter(abs(.std.resid) > 3) -> ofav_heat_outliers
ssid_model_heat_metrics %>% 
  filter(abs(.std.resid) > 3) -> ssid_heat_outliers

# remove outliers
# ipam_ofav_heat <- ipam_ofav_heat[!rownames(ipam_ofav_heat) %in% ofav_heat_outliers$.rownames,]
ipam_ssid_heat <- ipam_ssid_heat[!rownames(ipam_ssid_heat) %in% ssid_heat_outliers$.rownames,]

#### ANCOVA with heat ####

ofav_heat_aov <- ipam_ofav_heat %>% 
  anova_test(ratio ~ time + treatment*site)
get_anova_table(ofav_heat_aov) # all factors and treatment:site interaction significant
capture.output(get_anova_table(ofav_heat_aov), file = "urban ipam ofav heat ancova.txt")

ssid_heat_aov <- ipam_ssid_heat %>% 
  anova_test(ratio ~ time + treatment*site)
get_anova_table(ssid_heat_aov) # time, treatment, and site significant
capture.output(get_anova_table(ssid_heat_aov), file = "urban ipam ssid heat ancova.txt")


#### boxplots by species ####
fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a")

# creating dummy variables for stats labels on plots
ofav_control_time <- data.frame(time = 10,ratio = 1,lab = "time: F1,563=146.7, p<0.001*", pH = "low pH", temp="control temp")
ofav_control_treat <- data.frame(time = 10,ratio = 1,lab = "treat: F3,563=15.7, p<0.001*", pH = "low pH", temp="control temp")
ofav_control_site <- data.frame(time = 10,ratio = 1,lab = "site: F3,563=134.0, p<0.001*", pH = "low pH", temp="control temp")
ofav_control_int <- data.frame(time = 10,ratio = 1,lab = "int: F9,563=4.6, p<0.001*", pH = "low pH", temp="control temp")
ofav_heat_time <- data.frame(time = 10,ratio = 1,lab = "time: F1,594=348.6, p<0.001*", pH = "low pH", temp="high temp")
ofav_heat_treat <- data.frame(time = 10,ratio = 1,lab = "treat: F1,594=21.0, p<0.001*", pH = "low pH", temp="high temp")
ofav_heat_site <- data.frame(time = 10,ratio = 1,lab = "site: F3,594=98.7, p<0.001*", pH = "low pH", temp="high temp")
ofav_heat_int <- data.frame(time = 10,ratio = 1,lab = "int: F3,594=9.0, p<0.001*", pH = "low pH", temp="high temp")

ofav <- ggboxplot(ipam_ofav,
                  x = "time",
                  y = "ratio",
                  fill = "site",
                  palette=fill.color, 
                  width = 0.7,
                  size = 0.75,
                  title = "Ofav",
                  legend = "bottom") + 
  facet_grid(temp ~ pH) +
  geom_text(data = ofav_control_time,label = ofav_control_time$lab, aes(x=12.5, y=1)) +
  geom_text(data = ofav_control_treat,label = ofav_control_treat$lab, aes(x=12.5, y=0.925)) +
  geom_text(data = ofav_control_site,label = ofav_control_site$lab, aes(x=12.5, y=0.85)) +
  geom_text(data = ofav_control_int,label = ofav_control_int$lab, aes(x=12.5, y=0.775)) +
  geom_text(data = ofav_heat_time,label = ofav_heat_time$lab, aes(x=12.5, y=1)) +
  geom_text(data = ofav_heat_treat,label = ofav_heat_treat$lab, aes(x=12.5, y=0.925)) +
  geom_text(data = ofav_heat_site,label = ofav_heat_site$lab, aes(x=12.5, y=0.85)) +
  geom_text(data = ofav_heat_int,label = ofav_heat_int$lab, aes(x=12.5, y=0.775)) +
  xlab("Day") + 
  ylab("Proportional Fv/Fm") 
ofav

ggsave("urban ipam ofav.pdf", ofav, width=12, height=6,dpi = 300)


# creating dummy variables for stats labels on plots
ssid_control_time <- data.frame(time = 10,ratio = 1,lab = "time: F1,644=35.0, p<0.001*", pH = "low pH", temp="control temp")
ssid_control_treat <- data.frame(time = 10,ratio = 1,lab = "treat: F3,644=115.1, p<0.001*", pH = "low pH", temp="control temp")
ssid_control_site <- data.frame(time = 10,ratio = 1,lab = "site: F3,644=1.4, p=0.24", pH = "low pH", temp="control temp")
ssid_heat_time <- data.frame(time = 10,ratio = 1,lab = "time: F1,697=6.7, p=0.01*", pH = "low pH", temp="high temp")
ssid_heat_treat <- data.frame(time = 10,ratio = 1,lab = "treat: F1,697=76.1, p<0.001*", pH = "low pH", temp="high temp")
ssid_heat_site <- data.frame(time = 10,ratio = 1,lab = "site: F3,697=5.2, p=0.001*", pH = "low pH", temp="high temp")

ssid <- ggboxplot(ipam_ssid,
                  x = "time",
                  y = "ratio",
                  fill = "site",
                  palette=fill.color, 
                  width = 0.7,
                  size = 0.75,
                  title = "Ssid",
                  legend = "bottom") + 
  facet_grid(temp ~ pH) +
  geom_text(data = ssid_control_time,label = ssid_control_time$lab, aes(x=12.5, y=1)) +
  geom_text(data = ssid_control_treat,label = ssid_control_treat$lab, aes(x=12.5, y=0.9)) +
  geom_text(data = ssid_control_site,label = ssid_control_site$lab, aes(x=12.5, y=0.8)) +
  geom_text(data = ssid_heat_time,label = ssid_heat_time$lab, aes(x=12.5, y=1)) +
  geom_text(data = ssid_heat_treat,label = ssid_heat_treat$lab, aes(x=12.5, y=0.9)) +
  geom_text(data = ssid_heat_site,label = ssid_heat_site$lab, aes(x=12.5, y=0.8)) +
  xlab("Day") + 
  ylab("Proportional Fv/Fm")
ssid

ggsave("urban ipam ssid.pdf", ssid, width=12, height=6,dpi = 300)

