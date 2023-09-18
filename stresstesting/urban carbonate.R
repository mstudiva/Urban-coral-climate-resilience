#### packages ####

library(tidyverse)
library(seacarb) #install.packages("seacarb")
library(ggplot2)
library(ggpubr)
library(stringr)
library(rcompanion)


#### data import ####

tank <- read.csv("urban carbonate treatments.csv", head=T) # weekly tank samples for treatment means
tank$species <- as.factor(tank$species)
tank$treatment <- as.factor(tank$treatment)
tank$TA_mol <- tank$TA/1000000 # converting TA umol/kg to mol/kg
str(tank)

incubations <- read.csv("urban carbonate incubations.csv", head=T) # incubation samples for alkalinity anomaly/calcification
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

write.csv(tank_means, file = "urban carbonate treatment means.csv")


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

write.csv(incubations_netmeans, file = "urban carbonate incubations net means.csv")


#### incubation stats ####

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

pdf("urban incubations ofav normality.pdf") # histograms and Q-Q plots
par(mfrow=c(2,2))
hist(incubations_ofav_day$calcificationadj)
hist(incubations_ofav_night$calcificationadj)
qqnorm(incubations_ofav_day$calcificationadj)
qqline(incubations_ofav_day$calcificationadj)
qqnorm(incubations_ofav_night$calcificationadj)
qqline(incubations_ofav_night$calcificationadj)
dev.off()

pdf("urban incubations ssid normality.pdf") # histograms and Q-Q plots
par(mfrow=c(2,2))
hist(incubations_ssid_day$calcificationadj)
hist(incubations_ssid_night$calcificationadj)
qqnorm(incubations_ssid_day$calcificationadj)
qqline(incubations_ssid_day$calcificationadj)
qqnorm(incubations_ssid_night$calcificationadj)
qqline(incubations_ssid_night$calcificationadj)
dev.off()

anova_ofav_day <- aov(calcificationadj ~ treatment*site+site/genotype, data=incubations_ofav_day)
summary(anova_ofav_day) # no factors significant, genotype marginal
capture.output(summary(anova_ofav_day), file = "urban carbonate ofav day anova.txt")

anova_ofav_night <- aov(calcificationadj ~ treatment*site+site/genotype, data=incubations_ofav_night)
summary(anova_ofav_night) # treatment marginal
capture.output(summary(anova_ofav_night), file = "urban carbonate ofav night anova.txt")

anova_ssid_day <- aov(calcificationadj ~ treatment*site+site/genotype, data=incubations_ssid_day)
summary(anova_ssid_day) # treatment and site significant
capture.output(summary(anova_ssid_day), file = "urban carbonate ssid day anova.txt")

anova_ssid_night <- aov(calcificationadj ~ treatment*site+site/genotype, data=incubations_ssid_night)
summary(anova_ssid_night) # treatment and site significant
capture.output(summary(anova_ssid_night), file = "urban carbonate ssid night anova.txt")

# tukey_ofav_day <- TukeyHSD(anova_ofav_day) # pairwise comparisons
# tukey_ofav_day$`treatment:site`
# capture.output(tukey_ofav_day$`treatment:site`, file = "urban carbonate ofav day tukey.txt")

# tukey_ofav_night <- TukeyHSD(anova_ofav_night) # pairwise comparisons
# tukey_ofav_night$`treatment:site`
# capture.output(tukey_ofav_night$`treatment:site`, file = "urban carbonate ofav night tukey.txt")

tukey_ssid_day <- TukeyHSD(anova_ssid_day) # pairwise comparisons
tukey_ssid_day$`treatment:site`
capture.output(tukey_ssid_day$`treatment:site`, file = "urban carbonate ssid day tukey.txt")

tukey_ssid_night <- TukeyHSD(anova_ssid_night) # pairwise comparisons
tukey_ssid_night$`treatment:site`
capture.output(tukey_ssid_night$`treatment:site`, file = "urban carbonate ssid night tukey.txt")

# letters_ofav_day <- data.frame(tukey_ofav_day$`treatment:site`) # significance letters for plot
# letters_ofav_day$Var <- rownames(letters_ofav_day)
# names(letters_ofav_day)[5] <- "comparison"
# letters_ofav_day$comparison = str_replace_all(letters_ofav_day$comparison,":","_")
# letters_ofav_day$p.adj[is.na(letters_ofav_day$p.adj)] <- 1
# letters_ofav_day

# letters_ofav_night <- data.frame(tukey_ofav_night$`treatment:site`) # significance letters for plot
# letters_ofav_night$Var <- rownames(letters_ofav_night)
# names(letters_ofav_night)[5] <- "comparison"
# letters_ofav_night$comparison = str_replace_all(letters_ofav_night$comparison,":","_")
# letters_ofav_night$p.adj[is.na(letters_ofav_night$p.adj)] <- 1
# letters_ofav_night

letters_ssid_day <- data.frame(tukey_ssid_day$`treatment:site`) # significance letters for plot
letters_ssid_day$Var <- rownames(letters_ssid_day)
names(letters_ssid_day)[5] <- "comparison"
letters_ssid_day$comparison = str_replace_all(letters_ssid_day$comparison,":","_")
letters_ssid_day$p.adj[is.na(letters_ssid_day$p.adj)] <- 1
letters_ssid_day

letters_ssid_night <- data.frame(tukey_ssid_night$`treatment:site`) # significance letters for plot
letters_ssid_night$Var <- rownames(letters_ssid_night)
names(letters_ssid_night)[5] <- "comparison"
letters_ssid_night$comparison = str_replace_all(letters_ssid_night$comparison,":","_")
letters_ssid_night$p.adj[is.na(letters_ssid_night$p.adj)] <- 1
letters_ssid_night

# cld_ofav_day <- cldList(p.adj ~ comparison, data = letters_ofav_day, threshold = 0.05) # compact letter display for plot
# cld_ofav_day %>%
#   separate(Group, c("treatment", "site"), sep="_") %>%
#   mutate(across('treatment', str_replace, 'controlpH', 'control pH')) %>%
#   mutate(across('treatment', str_replace, 'lowpH', 'low pH')) %>%
#   mutate(across('site', str_replace, 'MacArthurNorth', 'MacArthur North')) %>%
#   mutate(across('site', str_replace, 'StarIsland', 'Star Island')) -> cld_ofav_day
# cld_ofav_day

# cld_ofav_night <- cldList(p.adj ~ comparison, data = letters_ofav_night, threshold = 0.05) # compact letter display for plot
# cld_ofav_night %>%
#   separate(Group, c("treatment", "site"), sep="_") %>%
#   mutate(across('treatment', str_replace, 'controlpH', 'control pH')) %>%
#   mutate(across('treatment', str_replace, 'lowpH', 'low pH')) %>%
#   mutate(across('site', str_replace, 'MacArthurNorth', 'MacArthur North')) %>%
#   mutate(across('site', str_replace, 'StarIsland', 'Star Island')) -> cld_ofav_night
# cld_ofav_night

cld_ssid_day <- cldList(p.adj ~ comparison, data = letters_ssid_day, threshold = 0.05) # compact letter display for plot
cld_ssid_day %>%
  separate(Group, c("treatment", "site"), sep="_") %>%
  mutate(across('treatment', str_replace, 'controlpH', 'control pH')) %>%
  mutate(across('treatment', str_replace, 'lowpH', 'low pH')) %>%
  mutate(across('site', str_replace, 'MacArthurNorth', 'MacArthur North')) %>%
  mutate(across('site', str_replace, 'StarIsland', 'Star Island')) -> cld_ssid_day
cld_ssid_day

cld_ssid_night <- cldList(p.adj ~ comparison, data = letters_ssid_night, threshold = 0.05) # compact letter display for plot
cld_ssid_night %>%
  separate(Group, c("treatment", "site"), sep="_") %>%
  mutate(across('treatment', str_replace, 'controlpH', 'control pH')) %>%
  mutate(across('treatment', str_replace, 'lowpH', 'low pH')) %>%
  mutate(across('site', str_replace, 'MacArthurNorth', 'MacArthur North')) %>%
  mutate(across('site', str_replace, 'StarIsland', 'Star Island')) -> cld_ssid_night
cld_ssid_night


#### incubation plots ####

fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a","grey") # custom color palette


stats_ofav_day_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,41 = 1.3, p = 0.3") # creating dummy variables for stats labels on plot
stats_ofav_day_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,41 = 0.8, p = 0.5")
stats_ofav_day_genotype <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "geno: F17,41 = 2.0, p = 0.08")

stats_ofav_night_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,39 = 3.5, p = 0.08") # creating dummy variables for stats labels on plot
stats_ofav_night_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,39 = 0.2, p = 0.9")
stats_ofav_night_genotype <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "geno: F17,39 = 0.7, p = 0.8")

stats_ssid_day_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,39 = 14.6, p = 0.002*") # creating dummy variables for stats labels on plot
stats_ssid_day_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,39 = 8.3, p = 0.001*")
stats_ssid_day_genotype <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "geno: F16,39 = 1.5, p = 0.2")

stats_ssid_night_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,37 = 29.6, p < 0.001*") # creating dummy variables for stats labels on plot
stats_ssid_night_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,37 = 4.9, p = 0.02*")
stats_ssid_night_genotype <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "geno: F16,37 = 1.9, p = 0.1")

ofav_day <- ggboxplot(incubations_ofav_day,
                  x = "site",
                  y = "calcificationadj",
                  fill = "site",
                  width = 0.7,
                  size = 0.75,
                  title = "Orbicella faveolata Day",
                  legend = "none") +
  facet_grid(~treatment) +
  xlab(element_blank()) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
  geom_text(data = stats_ofav_day_ph,label = stats_ofav_day_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ofav_day_site,label = stats_ofav_day_site$lab, aes(x=2, y=0.085)) +
  geom_text(data = stats_ofav_day_genotype,label = stats_ofav_day_genotype$lab, aes(x=2, y=0.08)) +
  # geom_text(data=cld_ofav_day, aes(x = site, y=-0.025, label=Letter)) +
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
  xlab(element_blank()) +
  ylab(element_blank()) +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(),axis.text.y = element_blank()) +
  geom_text(data = stats_ofav_night_ph,label = stats_ofav_night_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ofav_night_site,label = stats_ofav_night_site$lab, aes(x=2, y=0.085)) +
  geom_text(data = stats_ofav_night_genotype,label = stats_ofav_night_genotype$lab, aes(x=2, y=0.08)) +
  # geom_text(data=cld_ofav_night, aes(x = site, y=-0.025, label=Letter)) +
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
  xlab(element_blank()) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = stats_ssid_day_ph,label = stats_ssid_day_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ssid_day_site,label = stats_ssid_day_site$lab, aes(x=2, y=0.085)) +
  geom_text(data = stats_ssid_day_genotype,label = stats_ssid_day_genotype$lab, aes(x=2, y=0.08)) +
  geom_text(data=cld_ssid_day, aes(x = site, y=-0.025, label=Letter)) +
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
  xlab(element_blank()) +
  ylab(element_blank()) +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_blank()) +
  geom_text(data = stats_ssid_night_ph,label = stats_ssid_night_ph$lab, aes(x=2, y=0.09)) +
  geom_text(data = stats_ssid_night_site,label = stats_ssid_night_site$lab, aes(x=2, y=0.085)) +
  geom_text(data = stats_ssid_night_genotype,label = stats_ssid_night_genotype$lab, aes(x=2, y=0.08)) +
  geom_text(data=cld_ssid_night, aes(x = site, y=-0.025, label=Letter)) +
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

ggsave("urban carbonate incubations.pdf", incubations_species, width=16, height=8,dpi = 300)


#### incubation day/night means ####

incubations_net_carb %>% # generating mean adjusted calcification per sample (average day/night)
  group_by(sampleID, site, species,treatment,genotype)  %>% 
  summarise(calcificationavg = mean(calcificationadj)) -> incubations_meanCalc

incubations_meanCalc %>% # generating summary statistics for average adjusted calcification
  group_by(site, species,treatment)  %>% 
  summarise(meanCalc = mean(calcificationavg), sdCalc = sd(calcificationavg), nCalc = n(), seCalc = sdCalc/sqrt(nCalc)) -> incubations_mean_netCalc

write.csv(incubations_mean_netCalc, file = "urban carbonate incubations daynight calcification.csv")

incubations_mean_subset <- subset(incubations_meanCalc, site!="control") # filtering out skeletal controls

incubations_mean_ofav <- incubations_mean_subset[incubations_mean_subset$species == "Ofav",] # filtering dataframe by species 
incubations_mean_ssid <- incubations_mean_subset[incubations_mean_subset$species == "Ssid",] # filtering dataframe by species 

# identifying outliers
outliers_ofav <- boxplot(incubations_mean_ofav$calcificationavg, plot=T)$out # no outliers
outliers_ssid <- boxplot(incubations_mean_ssid$calcificationavg, plot=T)$out # no outliers

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(incubations_mean_ofav$calcificationavg) # normal
shapiro.test(incubations_mean_ssid$calcificationavg) # normal

pdf("urban incubations daynight normality.pdf") # histograms and Q-Q plots
par(mfrow=c(2,2))
hist(incubations_mean_ofav$calcificationavg)
hist(incubations_mean_ssid$calcificationavg)
qqnorm(incubations_mean_ofav$calcificationavg)
qqline(incubations_mean_ofav$calcificationavg)
qqnorm(incubations_mean_ssid$calcificationavg)
qqline(incubations_mean_ssid$calcificationavg)
dev.off()

anova_mean_ofav <- aov(calcificationavg ~ treatment*site+site/genotype, data=incubations_mean_ofav)
summary(anova_mean_ofav) # no factors significant
capture.output(summary(anova_mean_ofav), file = "urban carbonate mean ofav anova.txt")

anova_mean_ssid <- aov(calcificationavg ~ treatment*site+site/genotype, data=incubations_mean_ssid)
summary(anova_mean_ssid) # treatment and site significant, genotype marginal
capture.output(summary(anova_mean_ssid), file = "urban carbonate mean ssid anova.txt")

tukey_mean_ssid <- TukeyHSD(anova_mean_ssid) # pairwise comparisons
tukey_mean_ssid$`treatment:site`
capture.output(tukey_mean_ssid$`treatment:site`, file = "urban carbonate mean ssid tukey.txt")

letters_mean_ssid <- data.frame(tukey_mean_ssid$`treatment:site`) # significance letters for plot
letters_mean_ssid$Var <- rownames(letters_mean_ssid)
names(letters_mean_ssid)[5] <- "comparison"
letters_mean_ssid$comparison = str_replace_all(letters_mean_ssid$comparison,":","_")
letters_mean_ssid$p.adj[is.na(letters_mean_ssid$p.adj)] <- 1
letters_mean_ssid

cld_mean_ssid <- cldList(p.adj ~ comparison, data = letters_mean_ssid, threshold = 0.05) # compact letter display for plot
cld_mean_ssid %>%
  separate(Group, c("treatment", "site"), sep="_") %>%
  mutate(across('treatment', str_replace, 'controlpH', 'control pH')) %>%
  mutate(across('treatment', str_replace, 'lowpH', 'low pH')) %>%
  mutate(across('site', str_replace, 'MacArthurNorth', 'MacArthur North')) %>%
  mutate(across('site', str_replace, 'StarIsland', 'Star Island')) -> cld_mean_ssid
cld_mean_ssid

fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a","grey") # custom color palette

stats_mean_ofav_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,41 = 1.5, p = 0.2") # creating dummy variables for stats labels on plot
stats_mean_ofav_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,41 = 0.5, p = 0.7")
stats_mean_ofav_genotype <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "geno: F17,41 = 1.4, p = 0.3")

stats_mean_ssid_ph <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "pH: F1,39 = 17.3, p < 0.001*") # creating dummy variables for stats labels on plot
stats_mean_ssid_site <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "site: F3,39 = 8.6, p = 0.001*")
stats_mean_ssid_genotype <- data.frame(treatment = "control pH", growthadj = 0.6, lab = "geno: F16,39 = 2.1, p = 0.07")

ofav_daynight <- ggboxplot(incubations_mean_ofav,
                      x = "site",
                      y = "calcificationavg",
                      fill = "site",
                      width = 0.7,
                      size = 0.75,
                      title = "Orbicella faveolata Mean",
                      legend = "none") +
  facet_grid(~treatment) +
  xlab(element_blank()) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
  geom_text(data = stats_mean_ofav_ph,label = stats_mean_ofav_ph$lab, aes(x=2, y=0.08)) +
  geom_text(data = stats_mean_ofav_site,label = stats_mean_ofav_site$lab, aes(x=2, y=0.075)) +
  geom_text(data = stats_mean_ofav_genotype,label = stats_mean_ofav_genotype$lab, aes(x=2, y=0.07)) +
  # geom_text(data=cld_mean_ofav, aes(x = site, y=-0.025, label=Letter)) +
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
  xlab(element_blank()) +
  ylab("Calcification (mg cm-2 hr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = stats_mean_ssid_ph,label = stats_mean_ssid_ph$lab, aes(x=2, y=0.08)) +
  geom_text(data = stats_mean_ssid_site,label = stats_mean_ssid_site$lab, aes(x=2, y=0.075)) +
  geom_text(data = stats_mean_ssid_genotype,label = stats_mean_ssid_genotype$lab, aes(x=2, y=0.07)) +
  geom_text(data=cld_mean_ssid, aes(x = site, y=0.005, label=Letter)) +
  ylim(-0.01, 0.08)
ssid_daynight

incubations_daynight <-ggarrange(ofav_daynight, # both panels
                                 ssid_daynight,
                                heights = c(3.9,4),
                                widths = c(4),
                                ncol = 1,
                                nrow = 2)
incubations_daynight

ggsave("urban carbonate mean incubations.pdf", incubations_daynight, width=8, height=8,dpi = 300)
