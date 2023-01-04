#### packages ####

library(ggplot2)
library(ggpubr)
library(rcompanion)
library(MASS)
library(stringr)
library(survival)
library(survminer)
library(dplyr)


#### data import ####

fisher <- read.csv("fisher survivorship.csv", head=T)
str(fisher)
head(fisher)


#### survivorship by species ####

# create survival object using the Kaplan-Meier method
surv <- Surv(time = fisher$time, event = fisher$status)
surv

# run survival model for each species
fit_species <- survfit(surv ~ species, data = fisher)
summary(fit_species)

# Kaplan-Meier plots for each species
fill.color<-c("#43a2ca", "#a97c50","#addd8e", "#006837")

survival_species<-ggsurvplot(fit_species, data = fisher, pval = TRUE, xlab="Days", ylab="Survival Probability",
                                    conf.int = T, risk.table=T, palette=fill.color,
                                    break.time.by=31, xlim=c(0,372), risk.table.y.text = FALSE) + ggtitle("Species") 
survival_species

# hazard ratio by species
hazard_species <- coxph(surv ~ species, data = fisher)
summary(hazard_species)

# hazard ratio plot
hazplot_species <- ggforest(hazard_species, data = fisher)
hazplot_species


#### survivorship by size ####

# run survival model for both fragment size classes
fit_size <- survfit(surv ~ size, data = fisher)
summary(fit_size)

# Kaplan-Meier plots for each size
fill.color2<-c("grey75", "grey25")

survival_size<-ggsurvplot(fit_size, data = fisher, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color2,
                             break.time.by=31, xlim=c(0,372), risk.table.y.text = FALSE) + ggtitle("Size") 
survival_size

# hazard ratio by size
hazard_size <- coxph(surv ~ size, data = fisher)
summary(hazard_size)

# hazard ratio plot
hazplot_size <- ggforest(hazard_size, data = fisher)
hazplot_size


#### multiplot #####

# survivorship by species and size figure
fisher_multiplot<-ggarrange(survival_species$plot,
                                    survival_size$plot, 
                                    survival_species$table, 
                                    survival_size$table,
                                     hazplot_species,
                                    hazplot_size,
                                     heights = c(2, 0.5, 0.75),
                                     ncol = 2, nrow = 3)
fisher_multiplot

ggsave("fisher survivorship.pdf", fisher_multiplot, width=10, height=10,dpi = 300)
