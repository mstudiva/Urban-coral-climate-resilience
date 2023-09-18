#### packages ####

library(ggplot2)
library(ggpubr)
library(stringr)
library(rcompanion)
library(tidyr)
library(dplyr)


#### data import ####

bw <- read.csv(file="urban buoyant weight.csv", head=T)

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

# identifying outliers
outliers_ofav <- boxplot(bw_ofav$growthadj, plot=FALSE)$out # 1 outlier
bw_ofav <- bw_ofav[-which(bw_ofav$growthadj %in% outliers_ofav),]

outliers_ssid <- boxplot(bw_ssid$growthadj, plot=FALSE)$out # no outliers

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(bw_ofav$growthadj)
# normal

shapiro.test(bw_ssid$growthadj)
# not normal

pdf("urban bw normality.pdf")
par(mfrow=c(2,2))
hist(bw_ofav$growthadj)
hist(bw_ssid$growthadj)
qqnorm(bw_ofav$growthadj)
qqline(bw_ofav$growthadj)
qqnorm(bw_ssid$growthadj)
qqline(bw_ssid$growthadj)
dev.off()
# both species data appear normally distributed


#### statistical tests ####

# ANOVA
# anova <- aov(growthadj ~ pH*species*site+species/site/genotype, data=bw_subset)
# summary(anova)
# capture.output(summary(anova), file = "urban bw anova.txt")

anova_ofav <- aov(growthadj ~ pH*site+site/genotype, data=bw_ofav)
summary(anova_ofav)
# site and site:genotype significant
capture.output(summary(anova_ofav), file = "urban bw ofav anova.txt")

anova_ssid <- aov(growthadj ~ pH*site+site/genotype, data=bw_ssid)
summary(anova_ssid)
# site and pH significant
capture.output(summary(anova_ssid), file = "urban bw ssid anova.txt")

# Tukey post hoc tests
# tukey <- TukeyHSD(anova)
# tukey$`species:site`
# capture.output(tukey$`species:site`, file = "urban bw tukey.txt")

tukey_ofav <- TukeyHSD(anova_ofav)
tukey_ofav$`pH:site`
capture.output(tukey_ofav$`pH:site`, file = "urban bw ofav tukey.txt")

tukey_ssid <- TukeyHSD(anova_ssid)
tukey_ssid$`pH:site`
capture.output(tukey_ssid$`pH:site`, file = "urban bw ssid tukey.txt")

# creating dataframes of the pairwise comparisons needed for plots and doing a bit of table reformatting
# letters <- data.frame(tukey$`species:site`)
# letters$Var <- rownames(letters)
# names(letters)[5] <- "comparison"
# letters$comparison = str_replace_all(letters$comparison,":","_")
# letters$p.adj[is.na(letters$p.adj)] <- 1
# letters

letters_ofav <- data.frame(tukey_ofav$`pH:site`)
letters_ofav$Var <- rownames(letters_ofav)
names(letters_ofav)[5] <- "comparison"
letters_ofav$comparison = str_replace_all(letters_ofav$comparison,":","_")
letters_ofav$p.adj[is.na(letters_ofav$p.adj)] <- 1
letters_ofav

letters_ssid <- data.frame(tukey_ssid$`pH:site`)
letters_ssid$Var <- rownames(letters_ssid)
names(letters_ssid)[5] <- "comparison"
letters_ssid$comparison = str_replace_all(letters_ssid$comparison,":","_")
letters_ssid$p.adj[is.na(letters_ssid$p.adj)] <- 1
letters_ssid

# creates compact letter display of significant pairwise differences for figure
# cld <- cldList(p.adj ~ comparison, data = letters, threshold = 0.05)
# cld %>% separate(Group, c("species", "site")) -> cld
# cld

# subsetting for plotting
# cld_ofav <- cld[cld$species == "Ofav",] # filtering dataset by species
# cld_ssid <- cld[cld$species == "Ssid",]


cld_ofav <- cldList(p.adj ~ comparison, data = letters_ofav, threshold = 0.05)
cld_ofav %>%
  separate(Group, c("pH", "site"), sep="_") %>%
  mutate(across('pH', str_replace, 'controlpH', 'control pH')) %>%
  mutate(across('pH', str_replace, 'lowpH', 'low pH')) -> cld_ofav
cld_ofav

cld_ssid <- cldList(p.adj ~ comparison, data = letters_ssid, threshold = 0.05)
cld_ssid %>%
  separate(Group, c("pH", "site"), sep="_") %>%
  mutate(across('pH', str_replace, 'controlpH', 'control pH')) %>%
  mutate(across('pH', str_replace, 'lowpH', 'low pH')) -> cld_ssid
cld_ssid


#### plots ####

fill.color<-c("#018571","#80cdc1","#dfc27d","#a6611a","grey") # custom color palette

# creating dummy variables for stats labels on plot
# stats_site <- data.frame(species = "Ssid", growthadj = 0.6,lab = "site: F4,170=8.7, p<0.001*")
# stats_species <- data.frame(species = "Ssid", growthadj = 0.6,lab = "species: F1,170=14.1, p<0.001*")
# stats_genotype <- data.frame(species = "Ssid", growthadj = 0.6,lab = "genotype: F33,170=2.4, p<0.001*")
# stats_interaction <- data.frame(species = "Ssid", growthadj = 0.6,lab = "site:species: F4,170=8.3, p<0.001*")

stats_ofav_ph <- data.frame(pH = "control pH", growthadj = 0.6, lab = "pH: F1,82 = 0.2, p = 0.6")
stats_ofav_site <- data.frame(pH = "control pH", growthadj = 0.6, lab = "site: F3,82 = 18.2, p < 0.001*")
stats_ofav_genotype <- data.frame(pH = "control pH", growthadj = 0.6, lab = "geno: F17,82 = 4.1, p < 0.001*")

stats_ssid_ph <- data.frame(pH = "control pH", growthadj = 0.6, lab = "pH: F1,79 = 8.4, p = 0.005*")
stats_ssid_site <- data.frame(pH = "control pH", growthadj = 0.6, lab = "site: F3,79 = 4.3, p = 0.008*")
stats_ssid_genotype <- data.frame(pH = "control pH", growthadj = 0.6, lab = "geno: F16,79 = 1.6, p = 0.1")

# plot <- ggboxplot(bw_subset,
#                   x = "site",
#                   y = "growthadj",
#                   fill = "site",
#                   width = 0.7,
#                   size = 0.75,
#                   title = "Buoyant Weight",
#                   legend = "none") +
#   facet_grid(~species) +
#   # geom_text(data = stats_site,label = stats_site$lab, aes(x=4.55, y=-0.02)) +
# geom_text(data = stats_species,label = stats_species$lab, aes(x=4.35, y=-0.06)) +
# geom_text(data = stats_genotype,label = stats_genotype$lab, aes(x=4.29, y=-0.1)) +
# geom_text(data = stats_interaction,label = stats_interaction$lab, aes(x=4.25, y=-0.14)) +
#   xlab(element_blank()) +
#   ylab("growthadj Rate (g cm2 yr-1)") +
#   scale_fill_manual(values = fill.color) +
#   geom_hline(yintercept=0, linetype="dashed", color = "black") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   geom_text(data=cld_ofav, aes(x = site, y=-0.18, label=Letter)) +
#   geom_text(data=cld_ssid, aes(x = site, y=-0.18, label=Letter)) +
#   ylim(-0.18, 0.7)
# plot
# 
# ggsave("urban bw.pdf", plot, width=8, height=5,dpi = 300)

ofav <- ggboxplot(bw_ofav,
                  x = "site",
                  y = "growthadj",
                  fill = "site",
                  width = 0.7,
                  size = 0.75,
                  title = "Orbicella faveolata",
                  legend = "none") +
  facet_grid(~pH) +
  xlab(element_blank()) +
  ylab("Growth Rate (g cm-2 yr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank()) +
  geom_text(data = stats_ofav_ph,label = stats_ofav_ph$lab, aes(x=2, y=0.7)) +
  geom_text(data = stats_ofav_site,label = stats_ofav_site$lab, aes(x=2, y=0.65)) +
  geom_text(data = stats_ofav_genotype,label = stats_ofav_genotype$lab, aes(x=2, y=0.6)) +
  geom_text(data=cld_ofav, aes(x = site, y=-0.13, label=Letter)) +
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
  xlab(element_blank()) +
  ylab("Growth Rate (g cm-2 yr-1)") +
  scale_fill_manual(values = fill.color) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = stats_ssid_ph,label = stats_ssid_ph$lab, aes(x=2, y=0.7)) +
  geom_text(data = stats_ssid_site,label = stats_ssid_site$lab, aes(x=2, y=0.65)) +
  geom_text(data = stats_ssid_genotype,label = stats_ssid_genotype$lab, aes(x=2, y=0.6)) +
  geom_text(data=cld_ssid, aes(x = site, y=-0.13, label=Letter)) +
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

ggsave("urban bw.pdf", bw_species, width=8, height=8,dpi = 300)

