#### packages ####

library(tidyverse)
library(seacarb) #install.packages("seacarb")
library(ggplot2)
library(ggpubr)
library(stringr)
library(rcompanion)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)


#### data import ####

tank <- read.csv("../../data/calcification/urban carbonate treatments.csv", head=T) # weekly tank samples for treatment means
tank$species <- as.factor(tank$species)
tank$treatment <- as.factor(tank$treatment)
tank$TA_mol <- tank$TA/1000000 # converting TA umol/kg to mol/kg
str(tank)

incubations <- read.csv("../../data/calcification/urban carbonate incubations.csv", head=T) # incubation samples for alkalinity anomaly/calcification
incubations$incubation <- as.factor(incubations$incubation)
incubations$site <- factor(incubations$site, levels = c("Emerald", "Rainbow", "Star Island", "MacArthur North"))
incubations$genotype <- as.factor(incubations$genotype)
incubations$species <- as.factor(incubations$species)
incubations$treatment <- as.factor(incubations$treatment)
incubations$initialTA_mol <- incubations$initialTA/1000000 # converting TA umol/kg to mol/kg
incubations$finalTA_mol <- incubations$finalTA/1000000 # converting TA umol/kg to mol/kg
str(incubations)


#### treatment means ####

tank %>% # running the seacarb function to calculate other parameters from pH and TA
  mutate(carb(flag=8, pH, TA_mol, S=salinity, T=temp)) -> tank_carb

tank_carb %>% # generating summary statistics for TA
  filter(OA =="OA") %>%
  group_by(species,treatment)  %>% 
  summarise(meanTA = mean(TA), sdTA = sd(TA), nTA = n(), seTA = sdTA/sqrt(nTA)) -> tank_TA

tank_carb %>% # generating summary statistics for pH
  filter(OA =="OA") %>%
  group_by(species,treatment)  %>% 
  summarise(meanpH = mean(pH), sdpH = sd(pH), npH = n(), sepH = sdpH/sqrt(npH)) -> tank_pH

tank_carb %>% # generating summary statistics for DIC
  filter(OA =="OA") %>%
  group_by(species,treatment)  %>% 
  summarise(meanDIC = mean(DIC), sdDIC = sd(DIC), nDIC = n(), seDIC = sdDIC/sqrt(nDIC)) -> tank_DIC

tank_carb %>% # generating summary statistics for pCO2
  filter(OA =="OA") %>%
  group_by(species,treatment)  %>% 
  summarise(meanpCO2 = mean(pCO2), sdpCO2 = sd(pCO2), npCO2 = n(), sepCO2 = sdpCO2/sqrt(npCO2)) -> tank_pCO2

tank_carb %>% # generating summary statistics for aragonite saturation state (Omega Aragonite)
  filter(OA =="OA") %>%
  group_by(species,treatment)  %>% 
  summarise(meanAr = mean(OmegaAragonite), sdAr = sd(OmegaAragonite), nAr = n(), seAr = sdAr/sqrt(nAr)) -> tank_Ar

tank_TA %>% # joining all the summary stats into one table
  left_join(tank_pH, by=c('species'='species', 'treatment'='treatment')) %>%
  left_join(tank_DIC, by=c('species'='species', 'treatment'='treatment')) %>%
  left_join(tank_pCO2, by=c('species'='species', 'treatment'='treatment')) %>%
  left_join(tank_Ar, by=c('species'='species', 'treatment'='treatment')) -> tank_means

write.csv(tank_means, file = "../../outputs/calcification/urban carbonate treatment means.csv")


#### incubation means ####

incubations %>% # running the seacarb function to calculate other parameters from initial pH and TA
  mutate(carb(flag=8, initialpH, initialTA_mol, S=salinity, T=temp)) -> incubations_initial_carb

incubations %>% # running the seacarb function to calculate other parameters from final pH and TA
  mutate(carb(flag=8, finalpH, finalTA_mol, S=salinity, T=temp)) -> incubations_final_carb

incubations_initial_carb %>% # pulls the initial and final DIC, pCO2, and aragonite saturation columns from the respective dataframes
  select(1:32, 41, 48, 50) %>%
  rename('initialpCO2'='pCO2','initialDIC'='DIC','initialAr'='OmegaAragonite') %>%
  left_join(select(incubations_final_carb, 1, 41, 48, 50), by='bottle') %>%
  rename('finalpCO2'='pCO2','finalDIC'='DIC','finalAr'='OmegaAragonite') -> incubations_carb

incubations_carb %>% # calculating net pH, DIC, pCO2, and aragonite saturation
  mutate(netDIC = finalDIC - initialDIC) %>%
  mutate(netpCO2 = finalpCO2 - initialpCO2) %>%
  mutate(netAr = finalAr - initialAr) -> incubations_net_carb
  
incubations_net_carb %>% # generating summary statistics for net TA
  group_by(incubation, site, species, treatment)  %>% 
  summarise(meanTA = mean(netTA), sdTA = sd(netTA), nTA = n(), seTA = sdTA/sqrt(nTA)) -> incubations_netTA

incubations_net_carb %>% # generating summary statistics for net pH
  group_by(incubation, site, species,treatment)  %>% 
  summarise(meanpH = mean(netpH), sdpH = sd(netpH), npH = n(), sepH = sdpH/sqrt(npH)) -> incubations_netpH

incubations_net_carb %>% # generating summary statistics for net DIC
  group_by(incubation, site, species,treatment)  %>% 
  summarise(meanDIC = mean(netDIC), sdDIC = sd(netDIC), nDIC = n(), seDIC = sdDIC/sqrt(nDIC)) -> incubations_netDIC

incubations_net_carb %>% # generating summary statistics for net pCO2
  group_by(incubation, site, species,treatment)  %>% 
  summarise(meanpCO2 = mean(netpCO2), sdpCO2 = sd(netpCO2), npCO2 = n(), sepCO2 = sdpCO2/sqrt(npCO2)) -> incubations_netpCO2

incubations_net_carb %>% # generating summary statistics for net aragonite saturation state (Omega Aragonite)
  group_by(incubation, site, species,treatment)  %>% 
  summarise(meanAr = mean(netAr), sdAr = sd(netAr), nAr = n(), seAr = sdAr/sqrt(nAr)) -> incubations_netAr

incubations_net_carb %>% # generating summary statistics for net adjusted calcification
  group_by(incubation, site, species,treatment)  %>% 
  summarise(meanCalc = mean(calcificationadj), sdCalc = sd(calcificationadj), nCalc = n(), seCalc = sdCalc/sqrt(nCalc)) -> incubations_netCalc


incubations_netTA %>% # joining all the summary stats into one table
  left_join(incubations_netpH, by=c('incubation'='incubation','site'='site','species'='species', 'treatment'='treatment')) %>%
  left_join(incubations_netDIC, by=c('incubation'='incubation','site'='site','species'='species', 'treatment'='treatment')) %>%
  left_join(incubations_netpCO2, by=c('incubation'='incubation','site'='site','species'='species', 'treatment'='treatment')) %>%
  left_join(incubations_netAr, by=c('incubation'='incubation','site'='site','species'='species', 'treatment'='treatment')) %>%
  left_join(incubations_netCalc, by=c('incubation'='incubation','site'='site','species'='species', 'treatment'='treatment')) -> incubations_netmeans

write.csv(incubations_netmeans, file = "../../outputs/calcification/urban carbonate incubations net means.csv")


#### incubation normality and assumption testing ####

incubations_subset <- subset(incubations_net_carb, site!="control") # filtering out skeletal controls

incubations_ofav_day <- incubations_subset[incubations_subset$species == "Ofav" & incubations_subset$incubation == "day",] # filtering dataframe by species and incubation type
incubations_ofav_night <- incubations_subset[incubations_subset$species == "Ofav" & incubations_subset$incubation == "night",] 
incubations_ssid_day <- incubations_subset[incubations_subset$species == "Ssid" & incubations_subset$incubation == "day",]
incubations_ssid_night <- incubations_subset[incubations_subset$species == "Ssid" & incubations_subset$incubation == "night",]

# identifying outliers
outliers_ofav_day <- boxplot(incubations_ofav_day$calcificationadj, plot=T)$out # no outliers
outliers_ofav_night <- boxplot(incubations_ofav_night$calcificationadj, plot=T)$out # 2 outliers
incubations_ofav_night <- incubations_ofav_night[-which(incubations_ofav_night$calcificationadj %in% outliers_ofav_night),] # outlier removal

outliers_ssid_day <- boxplot(incubations_ssid_day$calcificationadj, plot=T)$out # no outliers
outliers_ssid_night <- boxplot(incubations_ssid_night$calcificationadj, plot=T)$out # 2 outliers
incubations_ssid_night <- incubations_ssid_night[-which(incubations_ssid_night$calcificationadj %in% outliers_ssid_night),]

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(incubations_ofav_day$calcificationadj) # normal
shapiro.test(incubations_ofav_night$calcificationadj) # normal

shapiro.test(incubations_ssid_day$calcificationadj) # normal
shapiro.test(incubations_ssid_night$calcificationadj) # normal

pdf("../../outputs/calcification/urban incubations ofav normality.pdf") # histograms and Q-Q plots
par(mfrow=c(2,2))
hist(incubations_ofav_day$calcificationadj)
hist(incubations_ofav_night$calcificationadj)
qqnorm(incubations_ofav_day$calcificationadj)
qqline(incubations_ofav_day$calcificationadj)
qqnorm(incubations_ofav_night$calcificationadj)
qqline(incubations_ofav_night$calcificationadj)
dev.off()

pdf("../../outputs/calcification/urban incubations ssid normality.pdf") # histograms and Q-Q plots
par(mfrow=c(2,2))
hist(incubations_ssid_day$calcificationadj)
hist(incubations_ssid_night$calcificationadj)
qqnorm(incubations_ssid_day$calcificationadj)
qqline(incubations_ssid_day$calcificationadj)
qqnorm(incubations_ssid_night$calcificationadj)
qqline(incubations_ssid_night$calcificationadj)
dev.off()


#### incubation linear mixed-effect models ####

# O. faveolata day
# full model, tank nested within pH, genotype nested within site
lm_ofav_day <- lmer(calcificationadj ~ treatment * site +
                  (1 | site:genotype),
                data = incubations_ofav_day, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ofav_day_noint <- lmer(calcificationadj ~ treatment + site +
                        (1 | site:genotype),
                      data = incubations_ofav_day, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ofav_day, REML = FALSE), update(lm_ofav_day_noint, REML = FALSE)) # p > 0.05 but AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ofav_day_noG <- lm(calcificationadj ~ treatment + site,
                    data = incubations_ofav_day)

# LRT of genotype effect
anova(update(lm_ofav_day_noint, REML = FALSE), update(lm_ofav_day_noG)) # p > 0.05 and AIC about the same, so keep genotype

# model outputs
summary(lm_ofav_day_noint)
anova(lm_ofav_day_noint) # no significant factors
VarCorr(lm_ofav_day_noint)

capture.output(anova(lm_ofav_day_noint), file = "../../outputs/calcification/urban incubations ofav day lme.txt")

# pairwise tests
# emm_ofav_day <- emmeans(lm_ofav_day_noint, ~ site | treatment)
# pairs(emm_ofav_day, adjust = "tukey")

# capture.output(pairs(emm_ofav_day, adjust = "tukey"), file = "../../outputs/calcification/urban incubations ofav day pairwise.txt")

# Create letters indicating significant differences for plot
# cld_ofav_day <- cld(emm_ofav_day, adjust = "tukey", Letters = letters, alpha = 0.05)
# cld_ofav_day


# O. faveolata night
# full model, tank nested within pH, genotype nested within site
lm_ofav_night <- lmer(calcificationadj ~ treatment * site +
                      (1 | site:genotype),
                    data = incubations_ofav_night, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ofav_night_noint <- lmer(calcificationadj ~ treatment + site +
                            (1 | site:genotype),
                          data = incubations_ofav_night, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ofav_night, REML = FALSE), update(lm_ofav_night_noint, REML = FALSE)) # p > 0.05 but AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ofav_night_noG <- lm(calcificationadj ~ treatment + site,
                      data = incubations_ofav_night)

# LRT of genotype effect
anova(update(lm_ofav_night_noint, REML = FALSE), update(lm_ofav_night_noG)) # p > 0.05 and AIC about the same, so keep genotype

# model outputs
summary(lm_ofav_night_noint)
anova(lm_ofav_night_noint) # treatment effect only 
VarCorr(lm_ofav_night_noint)

capture.output(anova(lm_ofav_night_noint), file = "../../outputs/calcification/urban incubations ofav night lme.txt")

# pairwise tests
emm_ofav_night <- emmeans(lm_ofav_night_noint, ~ treatment | site)
pairs(emm_ofav_night, adjust = "tukey")

capture.output(pairs(emm_ofav_night, adjust = "tukey"), file = "../../outputs/calcification/urban incubations ofav night pairwise.txt")

# Create letters indicating significant differences for plot
cld_ofav_night <- cld(emm_ofav_night, adjust = "tukey", Letters = letters, alpha = 0.05)
cld_ofav_night


# S. siderea day
# full model, tank nested within pH, genotype nested within site
lm_ssid_day <- lmer(calcificationadj ~ treatment * site +
                      (1 | site:genotype),
                    data = incubations_ssid_day, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ssid_day_noint <- lmer(calcificationadj ~ treatment + site +
                            (1 | site:genotype),
                          data = incubations_ssid_day, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ssid_day, REML = FALSE), update(lm_ssid_day_noint, REML = FALSE)) # p > 0.05 but AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ssid_day_noG <- lm(calcificationadj ~ treatment + site,
                      data = incubations_ssid_day)

# LRT of genotype effect
anova(update(lm_ssid_day_noint, REML = FALSE), update(lm_ssid_day_noG)) # p > 0.05 and AIC about the same, so keep genotype

# model outputs
summary(lm_ssid_day_noint)
anova(lm_ssid_day_noint) # treatment and site significant
VarCorr(lm_ssid_day_noint)

capture.output(anova(lm_ssid_day_noint), file = "../../outputs/calcification/urban incubations ssid day lme.txt")

# pairwise tests
emm_ssid_day <- emmeans(lm_ssid_day_noint, ~ site * treatment)
pairs(emm_ssid_day, adjust = "tukey")

capture.output(pairs(emm_ssid_day, adjust = "tukey"), file = "../../outputs/calcification/urban incubations ssid day pairwise.txt")

# Create letters indicating significant differences for plot
cld_ssid_day <- cld(emm_ssid_day, adjust = "tukey", Letters = letters, alpha = 0.05)
cld_ssid_day


# S. siderea night
# full model, tank nested within pH, genotype nested within site
lm_ssid_night <- lmer(calcificationadj ~ treatment * site +
                        (1 | site:genotype),
                      data = incubations_ssid_night, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ssid_night_noint <- lmer(calcificationadj ~ treatment + site +
                              (1 | site:genotype),
                            data = incubations_ssid_night, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ssid_night, REML = FALSE), update(lm_ssid_night_noint, REML = FALSE)) # p > 0.05 but AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ssid_night_noG <- lm(calcificationadj ~ treatment + site,
                        data = incubations_ssid_night)

# LRT of genotype effect
anova(update(lm_ssid_night_noint, REML = FALSE), update(lm_ssid_night_noG)) # p > 0.05 and AIC about the same, so keep genotype

# model outputs
summary(lm_ssid_night_noint)
anova(lm_ssid_night_noint) # treatment effect, marginally nonsignificant site effect 
VarCorr(lm_ssid_night_noint)

capture.output(anova(lm_ssid_night_noint), file = "../../outputs/calcification/urban incubations ssid night lme.txt")

# pairwise tests
emm_ssid_night <- emmeans(lm_ssid_night_noint, ~ treatment | site)
pairs(emm_ssid_night, adjust = "tukey")

capture.output(pairs(emm_ssid_night, adjust = "tukey"), file = "../../outputs/calcification/urban incubations ssid night pairwise.txt")

# Create letters indicating significant differences for plot
cld_ssid_night <- cld(emm_ssid_night, adjust = "tukey", Letters = letters, alpha = 0.05)
cld_ssid_night


#### incubation plots ####

fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a","grey") # custom color palette

stats_ofav_day_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,41 = 1.6, p = 0.2") # creating dummy variables for stats labels on plot
stats_ofav_day_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,41 = 0.4, p = 0.8")

stats_ofav_night_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,39 = 4.6, p = 0.04*") # creating dummy variables for stats labels on plot
stats_ofav_night_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,39 = 0.2, p = 0.9")

stats_ssid_day_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,39 = 13.8, p = 0.001*") # creating dummy variables for stats labels on plot
stats_ssid_day_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,39 = 5.6, p = 0.008*")

stats_ssid_night_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,37 = 28.6, p < 0.001*") # creating dummy variables for stats labels on plot
stats_ssid_night_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,37 = 2.7, p = 0.08")

ofav_day <- ggboxplot(incubations_ofav_day,
                  x = "site",
                  y = "calcificationadj",
                  fill = "site",
                  width = 0.7,
                  size = 0.75,
                  title = "Orbicella faveolata Day",
                  legend = "none") +
  facet_grid(~treatment) +
  xlab(NULL) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
  geom_text(data = stats_ofav_day_ph,label = stats_ofav_day_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ofav_day_site,label = stats_ofav_day_site$lab, aes(x=2, y=0.085)) +
  # geom_text(data=cld_ofav_day, aes(x = site, y=-0.025, label=.group)) +
  ylim(-0.025, 0.09)
ofav_day

ofav_night <- ggboxplot(incubations_ofav_night,
                      x = "site",
                      y = "calcificationadj",
                      fill = "site",
                      width = 0.7,
                      size = 0.75,
                      title = "Orbicella faveolata Night",
                      legend = "none") +
  facet_grid(~treatment) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),axis.text.y = element_blank()) +
  geom_text(data = stats_ofav_night_ph,label = stats_ofav_night_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ofav_night_site,label = stats_ofav_night_site$lab, aes(x=2, y=0.085)) +
  geom_text(data=cld_ofav_night, aes(x = site, y=-0.025, label=.group)) +
  ylim(-0.025, 0.09)
ofav_night

ssid_day <- ggboxplot(incubations_ssid_day,
                      x = "site",
                      y = "calcificationadj",
                      fill = "site",
                      width = 0.7,
                      size = 0.75,
                      title = "Siderastrea siderea Day",
                      legend = "none") +
  facet_grid(~treatment) +
  xlab(NULL) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = stats_ssid_day_ph,label = stats_ssid_day_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ssid_day_site,label = stats_ssid_day_site$lab, aes(x=2, y=0.085)) +
  geom_text(data=cld_ssid_day, aes(x = site, y=-0.025, label=.group)) +
  ylim(-0.025, 0.09)
ssid_day

ssid_night <- ggboxplot(incubations_ssid_night,
                        x = "site",
                        y = "calcificationadj",
                        fill = "site",
                        width = 0.7,
                        size = 0.75,
                        title = "Siderastrea siderea Night",
                        legend = "none") +
  facet_grid(~treatment) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_blank()) +
  geom_text(data = stats_ssid_night_ph,label = stats_ssid_night_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ssid_night_site,label = stats_ssid_night_site$lab, aes(x=2, y=0.085)) +
  geom_text(data=cld_ssid_night, aes(x = site, y=-0.025, label=.group)) +
  ylim(-0.025, 0.09)
ssid_night

incubations_species <-ggarrange(ofav_day, # both panels
                                ofav_night,
                                ssid_day,
                                ssid_night,
                      heights = c(4,4),
                      widths = c(4,3.5),
                      ncol = 2,
                      nrow = 2)
incubations_species

ggsave("../../outputs/calcification/urban carbonate incubations.pdf", incubations_species, width=16, height=8,dpi = 300)


#### incubation day/night means ####

incubations_net_carb %>% # generating mean adjusted calcification per sample (average day/night)
  group_by(sampleID, site, species,treatment,genotype)  %>% 
  summarise(calcificationavg = mean(calcificationadj)) -> incubations_meanCalc

write.csv(incubations_meanCalc, file = "../../outputs/calcification/urban carbonate incubations daynight calcification.csv")

incubations_meanCalc %>% # generating summary statistics for average adjusted calcification
  group_by(site, species,treatment)  %>% 
  summarise(meanCalc = mean(calcificationavg), sdCalc = sd(calcificationavg), nCalc = n(), seCalc = sdCalc/sqrt(nCalc)) -> incubations_mean_netCalc

write.csv(incubations_mean_netCalc, file = "../../outputs/calcification/urban carbonate incubations daynight calcification means.csv")

incubations_mean_subset <- subset(incubations_meanCalc, site!="control") # filtering out skeletal controls

incubations_mean_ofav <- incubations_mean_subset[incubations_mean_subset$species == "Ofav",] # filtering dataframe by species 
incubations_mean_ssid <- incubations_mean_subset[incubations_mean_subset$species == "Ssid",] # filtering dataframe by species 

# identifying outliers
outliers_ofav <- boxplot(incubations_mean_ofav$calcificationavg, plot=T)$out # no outliers
outliers_ssid <- boxplot(incubations_mean_ssid$calcificationavg, plot=T)$out # no outliers

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(incubations_mean_ofav$calcificationavg) # normal
shapiro.test(incubations_mean_ssid$calcificationavg) # normal

pdf("../../outputs/calcification/urban incubations daynight normality.pdf") # histograms and Q-Q plots
par(mfrow=c(2,2))
hist(incubations_mean_ofav$calcificationavg)
hist(incubations_mean_ssid$calcificationavg)
qqnorm(incubations_mean_ofav$calcificationavg)
qqline(incubations_mean_ofav$calcificationavg)
qqnorm(incubations_mean_ssid$calcificationavg)
qqline(incubations_mean_ssid$calcificationavg)
dev.off()


#### incubation day/night linear mixed-effect models ####

# O. faveolata 
# full model, tank nested within pH, genotype nested within site
lm_ofav <- lmer(calcificationavg ~ treatment * site +
                      (1 | site:genotype),
                    data = incubations_mean_ofav, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ofav_noint <- lmer(calcificationavg ~ treatment + site +
                            (1 | site:genotype),
                          data = incubations_mean_ofav, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ofav, REML = FALSE), update(lm_ofav_noint, REML = FALSE)) # p > 0.05 but AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ofav_noG <- lm(calcificationavg ~ treatment + site,
                      data = incubations_mean_ofav)

# LRT of genotype effect
anova(update(lm_ofav_noint, REML = FALSE), update(lm_ofav_noG)) # p > 0.05 and AIC about the same, so keep genotype

# model outputs
summary(lm_ofav_noint)
anova(lm_ofav_noint) # no significant factors
VarCorr(lm_ofav_noint)

capture.output(anova(lm_ofav_noint), file = "../../outputs/calcification/urban incubations ofav mean lme.txt")

# pairwise tests
# emm_ofav <- emmeans(lm_ofav_noint, ~ site * treatment)
# pairs(emm_ofav, adjust = "tukey")

# capture.output(pairs(emm_ofav, adjust = "tukey"), file = "../../outputs/calcification/urban incubations ofav mean pairwise.txt")

# Create letters indicating significant differences for plot
# cld_ofav <- cld(emm_ofav, adjust = "tukey", Letters = letters, alpha = 0.05)
# cld_ofav


# S. siderea night
# full model, tank nested within pH, genotype nested within site
lm_ssid <- lmer(calcificationavg ~ treatment * site +
                        (1 | site:genotype),
                      data = incubations_mean_ssid, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ssid_noint <- lmer(calcificationavg ~ treatment + site +
                              (1 | site:genotype),
                            data = incubations_mean_ssid, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ssid, REML = FALSE), update(lm_ssid_noint, REML = FALSE)) # p > 0.05 but AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ssid_noG <- lm(calcificationavg ~ treatment + site,
                        data = incubations_mean_ssid)

# LRT of genotype effect
anova(update(lm_ssid_noint, REML = FALSE), update(lm_ssid_noG)) # p > 0.05 and AIC about the same, so keep genotype

# model outputs
summary(lm_ssid_noint)
anova(lm_ssid_noint) # significant treatment and site effects 
VarCorr(lm_ssid_noint)

capture.output(anova(lm_ssid_noint), file = "../../outputs/calcification/urban incubations ssid mean lme.txt")

# pairwise tests
emm_ssid <- emmeans(lm_ssid_noint, ~ site * treatment)
pairs(emm_ssid, adjust = "tukey")

capture.output(pairs(emm_ssid, adjust = "tukey"), file = "../../outputs/calcification/urban incubations ssid mean pairwise.txt")

# Create letters indicating significant differences for plot
cld_ssid <- cld(emm_ssid, adjust = "tukey", Letters = letters, alpha = 0.05)
cld_ssid


#### incubation day/night plots ####

fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a","grey") # custom color palette

stats_mean_ofav_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,41 = 1.7, p = 0.2") # creating dummy variables for stats labels on plot
stats_mean_ofav_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,41 = 0.3, p = 0.8")

stats_mean_ssid_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,39 = 20.2, p < 0.001*") # creating dummy variables for stats labels on plot
stats_mean_ssid_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,39 = 4.1, p = 0.025*")

ofav_daynight <- ggboxplot(incubations_mean_ofav,
                      x = "site",
                      y = "calcificationavg",
                      fill = "site",
                      width = 0.7,
                      size = 0.75,
                      title = "Orbicella faveolata Mean",
                      legend = "none") +
  facet_grid(~treatment) +
  xlab(NULL) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
  geom_text(data = stats_mean_ofav_ph,label = stats_mean_ofav_ph$lab, aes(x=2, y=0.08)) +
  geom_text(data = stats_mean_ofav_site,label = stats_mean_ofav_site$lab, aes(x=2, y=0.075)) +
  # geom_text(data=cld_ofav, aes(x = site, y=-0.025, label=.group)) +
  ylim(-0.01, 0.08)
ofav_daynight

ssid_daynight <- ggboxplot(incubations_mean_ssid,
                           x = "site",
                           y = "calcificationavg",
                           fill = "site",
                           width = 0.7,
                           size = 0.75,
                           title = "Siderastrea siderea Mean",
                           legend = "none") +
  facet_grid(~treatment) +
  xlab(NULL) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = stats_mean_ssid_ph,label = stats_mean_ssid_ph$lab, aes(x=2, y=0.08)) +
  geom_text(data = stats_mean_ssid_site,label = stats_mean_ssid_site$lab, aes(x=2, y=0.075)) +
  geom_text(data=cld_ssid, aes(x = site, y=0.005, label=.group)) +
  ylim(-0.01, 0.08)
ssid_daynight

incubations_daynight <-ggarrange(ofav_daynight, # both panels
                                 ssid_daynight,
                                heights = c(3.9,4),
                                widths = c(4),
                                ncol = 1,
                                nrow = 2)
incubations_daynight

ggsave("../../outputs/calcification/urban carbonate mean incubations.pdf", incubations_daynight, width=8, height=8,dpi = 300)
