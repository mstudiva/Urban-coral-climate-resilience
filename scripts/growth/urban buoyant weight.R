#### packages ####

library(ggplot2)
library(ggpubr)
library(stringr)
library(rcompanion)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)


#### data import ####

bw <- read.csv(file="../../data/growth/urban buoyant weight.csv", head=T)

bw$tank <- as.factor(bw$tank)
bw$pH <- as.factor(bw$pH)
bw$site=factor(bw$site, levels=c("Emerald", "Rainbow", "Star", "MacN", "control")) 
bw$genotype <- as.factor(bw$genotype)

str(bw)
head(bw)


#### data normality and assumption testing ####

bw_subset <- subset(bw, site!="control")

bw_ofav <- bw_subset[bw_subset$species == "Ofav",] # filtering dataset by species
bw_ssid <- bw_subset[bw_subset$species == "Ssid",]
str(bw_ofav)
str(bw_ssid)

# identifying outliers with IQR > 1.5
outliers_ofav <- boxplot(bw_ofav$growthadj, plot=FALSE)$out # 1 outlier
bw_ofav <- bw_ofav[-which(bw_ofav$growthadj %in% outliers_ofav),]

outliers_ssid <- boxplot(bw_ssid$growthadj, plot=FALSE)$out # no outliers

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(bw_ofav$growthadj)
# normal

shapiro.test(bw_ssid$growthadj)
# not normal

pdf("../../outputs/growth/urban bw normality.pdf")
par(mfrow=c(2,2))
hist(bw_ofav$growthadj)
hist(bw_ssid$growthadj)
qqnorm(bw_ofav$growthadj)
qqline(bw_ofav$growthadj)
qqnorm(bw_ssid$growthadj)
qqline(bw_ssid$growthadj)
dev.off()
# both species data appear normally distributed


#### linear mixed-effect models ####

# O. faveolata
# full model, tank nested within pH, genotype nested within site
lm_ofav <- lmer(growthadj ~ pH * site + 
                  (1 | pH:tank) + 
                  (1 | site:genotype),
                data = bw_ofav, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ofav_noint <- lmer(growthadj ~ pH + site + 
                        (1 | pH:tank) + 
                        (1 | site:genotype),
                      data = bw_ofav, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ofav, REML = FALSE), update(lm_ofav_noint, REML = FALSE)) # p > 0.05 and AIC lower for noint model, so use reduced model

# reduced model, dropped genotype as factor
lm_ofav_noG <- lmer(growthadj ~ pH + site + 
                      (1 | pH:tank),
                    data = bw_ofav, REML = TRUE)

# LRT of genotype effect
anova(update(lm_ofav_noint, REML = FALSE), update(lm_ofav_noG, REML = FALSE)) # p < 0.05 and AIC much lower for noint model, so keep genotype

# sanity check: is there a pH effect when looking at growth just across tanks?
tank_means <- aggregate(growthadj ~ pH + tank, bw_ofav, mean)
summary(lm(growthadj ~ pH, data = tank_means)) # p > 0.05, no tank effect

# model outputs
summary(lm_ofav_noint)
anova(lm_ofav_noint) # site effect only 
VarCorr(lm_ofav_noint)

capture.output(anova(lm_ofav_noint), file = "../../outputs/growth/urban bw ofav lme.txt")

# pairwise site tests
emm_ofav <- emmeans(lm_ofav_noint, ~ site | pH)
pairs(emm_ofav, adjust = "tukey")

capture.output(pairs(emm_ofav, adjust = "tukey"), file = "../../outputs/growth/urban bw ofav pairwise.txt")

# Create letters indicating significant differences for plot
cld_ofav <- cld(emm_ofav, adjust = "tukey", Letters = letters, alpha = 0.05)
cld_ofav


# S. siderea
# full model, tank nested within pH, genotype nested within site
lm_ssid <- lmer(growthadj ~ pH * site + 
                  (1 | pH:tank) + 
                  (1 | site:genotype),
                data = bw_ssid, REML = TRUE)

# reduced model, no interaction of pH and site
lm_ssid_noint <- lmer(growthadj ~ pH + site + 
                        (1 | pH:tank) + 
                        (1 | site:genotype),
                      data = bw_ssid, REML = TRUE)

# likelihood ratio test of interactive vs non-interactive model
anova(update(lm_ssid, REML = FALSE), update(lm_ssid_noint, REML = FALSE)) # p = 0.06 and only marginally AIC lower for full model, so use reduced model 

# reduced model, dropped genotype as factor
lm_ssid_noG <- lmer(growthadj ~ pH + site + 
                      (1 | pH:tank),
                    data = bw_ssid, REML = TRUE)

# LRT of genotype effect
anova(update(lm_ssid_noint, REML = FALSE), update(lm_ssid_noG, REML = FALSE)) # p = 0.05 and AIC lower for noint model, so keep genotype

# sanity check: is there a pH effect when looking at growth just across tanks?
tank_means <- aggregate(growthadj ~ pH + tank, bw_ssid, mean)
summary(lm(growthadj ~ pH, data = tank_means)) # p > 0.05, no tank effect

# model outputs
summary(lm_ssid_noint)
anova(lm_ssid_noint) # marginally nonsignificant site effect
VarCorr(lm_ssid_noint)

capture.output(anova(lm_ssid_noint), file = "../../outputs/growth/urban bw ssid lme.txt")

# # pairwise site tests
# emm_ssid <- emmeans(lm_ssid_noint, ~ site | pH)
# pairs(emm_ssid, adjust = "tukey")
# 
# # Create letters indicating significant differences for plot
# cld_ssid <- cld(emm_ssid, adjust = "tukey", Letters = letters, alpha = 0.05)
# cld_ssid


#### plots ####

fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a","grey") # custom color palette

# creating dummy variables for stats labels on plot
stats_ofav_ph <- data.frame(pH = "control pH", growthadj = 0.6, lab = "pH: F1,82 = 0.02, p = 0.9")
stats_ofav_site <- data.frame(pH = "control pH", growthadj = 0.6, lab = "site: F3,82 = 4.3, p = 0.02*")

stats_ssid_ph <- data.frame(pH = "control pH", growthadj = 0.6, lab = "pH: F1,79 = 0.7, p = 0.5")
stats_ssid_site <- data.frame(pH = "control pH", growthadj = 0.6, lab = "site: F3,79 = 2.8, p = 0.07")

ofav <- ggboxplot(bw_ofav,
                  x = "site",
                  y = "growthadj",
                  fill = "site",
                  width = 0.7,
                  size = 0.75,
                  title = "Orbicella faveolata",
                  legend = "none") +
  facet_grid(~pH) +
  xlab(NULL) +
  ylab("Growth Rate (g cm-2 yr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
  geom_text(data = stats_ofav_ph,label = stats_ofav_ph$lab, aes(x=2, y=0.7)) +
  geom_text(data = stats_ofav_site,label = stats_ofav_site$lab, aes(x=2, y=0.65)) +
  geom_text(data=cld_ofav, aes(x = site, y=-0.13, label=.group)) +
  ylim(-0.13, 0.7)
ofav

ssid <- ggboxplot(bw_ssid,
                  x = "site",
                  y = "growthadj",
                  fill = "site",
                  width = 0.7,
                  size = 0.75,
                  title = "Siderastrea siderea",
                  legend = "none") +
  facet_grid(~pH) +
  xlab(NULL) +
  ylab("Growth Rate (g cm-2 yr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = stats_ssid_ph,label = stats_ssid_ph$lab, aes(x=2, y=0.7)) +
  geom_text(data = stats_ssid_site,label = stats_ssid_site$lab, aes(x=2, y=0.65)) +
  # geom_text(data=cld_ssid, aes(x = site, y=-0.13, label=.group)) +
  ylim(-0.13, 0.7)
ssid

# both panels
bw_species<-ggarrange(ofav,
                      ssid,
                      heights = c(3.9,4),
                      widths = c(4,4),
                      ncol = 1,
                      nrow = 2)
bw_species

ggsave("../../outputs/growth/urban bw.pdf", bw_species, width=8, height=8,dpi = 300)

