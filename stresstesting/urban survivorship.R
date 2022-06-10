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

urban <- read.csv("urban survivorship.csv", head=T)
str(urban)
head(urban)


#### survivorship by treatment ####

# subset dataframe by species
urban_ofav <- subset(urban, species=="Ofav")
urban_ofav$treat=factor(urban_ofav$treat, levels=c("CC","CH","LC","LH")) 
urban_ofav$site=factor(urban_ofav$site, levels=c("Emerald","Rainbow","Star","MacN")) 

urban_ssid <- subset(urban, species=="Ssid")
urban_ssid$treat=factor(urban_ssid$treat, levels=c("CC","CH","LC","LH")) 
urban_ssid$site=factor(urban_ssid$site, levels=c("Emerald","Rainbow","Star","MacN")) 

# create survival objects for each species (using the Kaplan-Meier method)
surv_ofav_treat <- Surv(time = urban_ofav$days, event = urban_ofav$status)
surv_ofav_treat

surv_ssid_treat <- Surv(time = urban_ssid$days, event = urban_ssid$status)
surv_ssid_treat

# run survival model for each species
fit_ofav_treat <- survfit(surv_ofav_treat ~ treat, data = urban_ofav)
summary(fit_ofav_treat)

fit_ssid_treat <- survfit(surv_ssid_treat ~ treat, data = urban_ssid)
summary(fit_ssid_treat)

# Kaplan-Meier plots for each species
fill.color<-c("#92c5de","#f4a582","#0571b0","#ca0020")

survival_ofav_treat<-ggsurvplot(fit_ofav_treat, data = urban_ofav, pval = FALSE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color, 
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata") 
survival_ofav_treat

survival_ssid_treat<-ggsurvplot(fit_ssid_treat, data = urban_ssid, pval = FALSE, xlab="Days", ylab="Survival Probability",
                          conf.int = T, risk.table=T, palette=fill.color, 
                          break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Siderastrea siderea") 
survival_ssid_treat

# survivorship by treatment figure (no stats)
urban_treatment_nostats_multiplot<-ggarrange(survival_ofav_treat$plot,
                                              survival_ssid_treat$plot, 
                                     survival_ofav_treat$table, 
                                     survival_ssid_treat$table,
                                     heights = c(2, 0.5),
                                     ncol = 2, nrow = 2)
urban_treatment_nostats_multiplot

ggsave("urban survivorship treatment nostats.pdf", urban_treatment_nostats_multiplot, width=10, height=8,dpi = 300)


#### survivorship by treatment subset ####

# Since mortality was not expected for control temp treatments, need to remove them for hazard ratio stats
urban_ofav_treat <- subset(urban_ofav, treat!="CC")
urban_ofav_treat <- subset(urban_ofav_treat, treat!="LC")
urban_ofav_treat$treat=factor(urban_ofav_treat$treat, levels=c("CH","LH")) 
urban_ofav_treat$site=factor(urban_ofav_treat$site, levels=c("Emerald","Rainbow","Star","MacN")) 

urban_ssid_treat <- subset(urban_ssid, treat!="CC")
urban_ssid_treat <- subset(urban_ssid_treat, treat!="LC")
urban_ssid_treat$treat=factor(urban_ssid_treat$treat, levels=c("CH","LH")) 
urban_ssid_treat$site=factor(urban_ssid_treat$site, levels=c("Emerald","Rainbow","Star","MacN")) 

# create survival objects for each species (using the Kaplan-Meier method)
surv_ofav_treat_sub <- Surv(time = urban_ofav_treat$days, event = urban_ofav_treat$status)
surv_ofav_treat_sub

surv_ssid_treat_sub <- Surv(time = urban_ssid_treat$days, event = urban_ssid_treat$status)
surv_ssid_treat_sub

# run survival model for each species
fit_ofav_treat_sub <- survfit(surv_ofav_treat_sub ~ treat, data = urban_ofav_treat)
summary(fit_ofav_treat_sub)

fit_ssid_treat_sub <- survfit(surv_ssid_treat_sub ~ treat, data = urban_ssid_treat)
summary(fit_ssid_treat_sub)

# Kaplan-Meier plots for each species
fill.color2<-c("#f4a582","#ca0020")

survival_ofav_treat_sub<-ggsurvplot(fit_ofav_treat_sub, data = urban_ofav_treat, pval = TRUE, xlab="Days", ylab="Survival Probability",
                                conf.int = T, risk.table=T, palette=fill.color2,
                                break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata") 
survival_ofav_treat_sub

survival_ssid_treat_sub<-ggsurvplot(fit_ssid_treat_sub, data = urban_ssid_treat, pval = TRUE, xlab="Days", ylab="Survival Probability",
                                conf.int = T, risk.table=T, palette=fill.color2,
                                break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Siderastrea siderea") 
survival_ssid_treat_sub

# hazard ratio by treatment
hazard_ofav_treat <- coxph(surv_ofav_treat_sub ~ treat, data = urban_ofav_treat)
summary(hazard_ofav_treat)

hazard_ssid_treat <- coxph(surv_ssid_treat_sub ~ treat, data = urban_ssid_treat)
summary(hazard_ssid_treat)

# hazard ratio plots
hazplot_ofav_treat <- ggforest(hazard_ofav_treat, data = urban_ofav_treat)
hazplot_ofav_treat

hazplot_ssid_treat <- ggforest(hazard_ssid_treat, data = urban_ssid_treat)
hazplot_ssid_treat

# survivorship by treatment figure
urban_treatment_multiplot<-ggarrange(survival_ofav_treat_sub$plot,
                                    survival_ssid_treat_sub$plot, 
                                 survival_ofav_treat_sub$table, 
                                 survival_ssid_treat_sub$table,
                                 hazplot_ofav_treat,
                                 hazplot_ssid_treat,
                                 heights = c(2, 0.5, 0.5),
                                 ncol = 2, nrow = 3)
urban_treatment_multiplot

ggsave("urban survivorship treatment.pdf", urban_treatment_multiplot, width=10, height=10,dpi = 300)


#### survivorship by site ####

# subsetting dataframes to separate out thermal stress treatments for each species
urban_ofav_ch <- subset(urban_ofav, treat=="CH")
urban_ofav_ch$site=factor(urban_ofav_ch$site, levels=c("Emerald","Rainbow","Star","MacN")) 

# Since sample sizes for this group are small, will combine sites to offshore/urban
urban_ofav_lh <- subset(urban_ofav, treat=="LH")
urban_ofav_lh <- urban_ofav_lh %>%
  mutate(type = case_when(
    site=="Emerald" ~ "Offshore",
    site=="Rainbow" ~ "Offshore",
    site=="Star" ~ "Urban",
    site=="MacN" ~ "Urban",
  ))
urban_ofav_lh$type=factor(urban_ofav_lh$type, levels=c("Offshore","Urban")) 

urban_ssid_ch <- subset(urban_ssid, treat=="CH")
urban_ssid_ch$site=factor(urban_ssid_ch$site, levels=c("Emerald","Rainbow","Star","MacN")) 

urban_ssid_lh <- subset(urban_ssid, treat=="LH")
urban_ssid_lh$site=factor(urban_ssid_lh$site, levels=c("Emerald","Rainbow","Star","MacN")) 

# create survival objects for each species/treatment (using the Kaplan-Meier method)
surv_ofav_ch <- Surv(time = urban_ofav_ch$days, event = urban_ofav_ch$status)
surv_ofav_ch

surv_ofav_lh <- Surv(time = urban_ofav_lh$days, event = urban_ofav_lh$status)
surv_ofav_lh

surv_ssid_ch <- Surv(time = urban_ssid_ch$days, event = urban_ssid_ch$status)
surv_ssid_ch

surv_ssid_lh <- Surv(time = urban_ssid_lh$days, event = urban_ssid_lh$status)
surv_ssid_lh

# run survival model for each species
fit_ofav_ch <- survfit(surv_ofav_ch ~ site, data = urban_ofav_ch)
summary(fit_ofav_ch)

fit_ofav_lh <- survfit(surv_ofav_lh ~ type, data = urban_ofav_lh)
summary(fit_ofav_lh)

fit_ssid_ch <- survfit(surv_ssid_ch ~ site, data = urban_ssid_ch)
summary(fit_ssid_ch)

fit_ssid_lh <- survfit(surv_ssid_lh ~ site, data = urban_ssid_lh)
summary(fit_ssid_lh)

# Kaplan-Meier plots for each species
fill.color3<-c("#018571","#80cdc1","#dfc27d","#a6611a")
fill.color4<-c("#018571","#a6611a")

survival_ofav_ch<-ggsurvplot(fit_ofav_ch, data = urban_ofav_ch, pval = TRUE, xlab="Days", ylab="Survival Probability",
                                    conf.int = T, risk.table=T, palette=fill.color3,
                                    break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata CH") 
survival_ofav_ch

survival_ofav_lh<-ggsurvplot(fit_ofav_lh, data = urban_ofav_lh, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color4,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata LH") 
survival_ofav_lh

survival_ssid_ch<-ggsurvplot(fit_ssid_ch, data = urban_ssid_ch, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color3,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Siderastrea siderea CH") 
survival_ssid_ch

survival_ssid_lh<-ggsurvplot(fit_ssid_lh, data = urban_ssid_lh, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color3,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Siderastrea siderea LH") 
survival_ssid_lh

# hazard ratio by treatment
hazard_ofav_ch <- coxph(surv_ofav_ch ~ site, data = urban_ofav_ch)
summary(hazard_ofav_ch)

hazard_ofav_lh <- coxph(surv_ofav_lh ~ type, data = urban_ofav_lh)
summary(hazard_ofav_lh)

hazard_ssid_ch <- coxph(surv_ssid_ch ~ site, data = urban_ssid_ch)
summary(hazard_ssid_ch)

hazard_ssid_lh <- coxph(surv_ssid_lh ~ site, data = urban_ssid_lh)
summary(hazard_ssid_lh)

# hazard ratio plots
hazplot_ofav_ch <- ggforest(hazard_ofav_ch, data = urban_ofav_ch)
hazplot_ofav_ch

hazplot_ofav_lh <- ggforest(hazard_ofav_lh, data = urban_ofav_lh) + coord_flip(ylim=c(0,1))
hazplot_ofav_lh

hazplot_ssid_ch <- ggforest(hazard_ssid_ch, data = urban_ssid_ch)
hazplot_ssid_ch

hazplot_ssid_lh <- ggforest(hazard_ssid_lh, data = urban_ssid_lh)
hazplot_ssid_lh

# survivorship by treatment/site figure
urban_ofav_site_multiplot<-ggarrange(survival_ofav_ch$plot,
                                     survival_ofav_lh$plot, 
                                     survival_ofav_ch$table, 
                                     survival_ofav_lh$table,
                                     hazplot_ofav_ch,
                                     hazplot_ofav_lh,
                                    heights = c(2, 0.5, 0.5),
                                    ncol = 2, nrow = 3)
urban_ofav_site_multiplot

ggsave("urban survivorship ofav site.pdf", urban_ofav_site_multiplot, width=10, height=10,dpi = 300)

urban_ssid_site_multiplot<-ggarrange(survival_ssid_ch$plot,
                                     survival_ssid_lh$plot, 
                                     survival_ssid_ch$table, 
                                     survival_ssid_lh$table,
                                     hazplot_ssid_ch,
                                     hazplot_ssid_lh,
                                     heights = c(2, 0.5, 0.5),
                                     ncol = 2, nrow = 3)
urban_ssid_site_multiplot

ggsave("urban survivorship ssid site.pdf", urban_ssid_site_multiplot, width=10, height=10,dpi = 300)


#### surviviorship by site type ####

# Since all site comparisons were insignificant, collapsing across sites to urban vs offshore
urban_ofav_site_ch <- urban_ofav_ch %>%
  mutate(type = case_when(
    site=="Emerald" ~ "Offshore",
    site=="Rainbow" ~ "Offshore",
    site=="Star" ~ "Urban",
    site=="MacN" ~ "Urban",
  ))
urban_ofav_site_ch$type=factor(urban_ofav_site_ch$type, levels=c("Offshore","Urban")) 

urban_ofav_site_lh <- urban_ofav_lh %>%
  mutate(type = case_when(
    site=="Emerald" ~ "Offshore",
    site=="Rainbow" ~ "Offshore",
    site=="Star" ~ "Urban",
    site=="MacN" ~ "Urban",
  ))
urban_ofav_site_lh$type=factor(urban_ofav_site_lh$type, levels=c("Offshore","Urban")) 

urban_ssid_site_ch <- urban_ssid_ch %>%
  mutate(type = case_when(
    site=="Emerald" ~ "Offshore",
    site=="Rainbow" ~ "Offshore",
    site=="Star" ~ "Urban",
    site=="MacN" ~ "Urban",
  ))
urban_ssid_site_ch$type=factor(urban_ssid_site_ch$type, levels=c("Offshore","Urban")) 

urban_ssid_site_lh <- urban_ssid_lh %>%
  mutate(type = case_when(
    site=="Emerald" ~ "Offshore",
    site=="Rainbow" ~ "Offshore",
    site=="Star" ~ "Urban",
    site=="MacN" ~ "Urban",
  ))
urban_ssid_site_lh$type=factor(urban_ssid_site_lh$type, levels=c("Offshore","Urban")) 

# create survival objects for each species/treatment (using the Kaplan-Meier method)
surv_ofav_site_ch <- Surv(time = urban_ofav_site_ch$days, event = urban_ofav_site_ch$status)
surv_ofav_site_ch

surv_ofav_site_lh <- Surv(time = urban_ofav_site_lh$days, event = urban_ofav_site_lh$status)
surv_ofav_site_lh

surv_ssid_site_ch <- Surv(time = urban_ssid_site_ch$days, event = urban_ssid_site_ch$status)
surv_ssid_site_ch

surv_ssid_site_lh <- Surv(time = urban_ssid_site_lh$days, event = urban_ssid_site_lh$status)
surv_ssid_site_lh

# run survival model for each species
fit_ofav_site_ch <- survfit(surv_ofav_site_ch ~ type, data = urban_ofav_site_ch)
summary(fit_ofav_site_ch)

fit_ofav_site_lh <- survfit(surv_ofav_site_lh ~ type, data = urban_ofav_site_lh)
summary(fit_ofav_site_lh)

fit_ssid_site_ch <- survfit(surv_ssid_site_ch ~ type, data = urban_ssid_site_ch)
summary(fit_ssid_site_ch)

fit_ssid_site_lh <- survfit(surv_ssid_site_lh ~ type, data = urban_ssid_site_lh)
summary(fit_ssid_site_lh)

# Kaplan-Meier plots for each species
survival_ofav_site_ch<-ggsurvplot(fit_ofav_site_ch, data = urban_ofav_site_ch, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color4,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata CH") 
survival_ofav_site_ch

survival_ofav_site_lh<-ggsurvplot(fit_ofav_site_lh, data = urban_ofav_site_lh, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color4,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Orbicella faveolata LH") 
survival_ofav_site_lh

survival_ssid_site_ch<-ggsurvplot(fit_ssid_site_ch, data = urban_ssid_site_ch, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color4,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Siderastrea siderea CH") 
survival_ssid_site_ch

survival_ssid_site_lh<-ggsurvplot(fit_ssid_site_lh, data = urban_ssid_site_lh, pval = TRUE, xlab="Days", ylab="Survival Probability",
                             conf.int = T, risk.table=T, palette=fill.color4,
                             break.time.by=5, xlim=c(0,36), risk.table.y.text = FALSE) + ggtitle("Siderastrea siderea LH") 
survival_ssid_site_lh

# hazard ratio by treatment
hazard_ofav_site_ch <- coxph(surv_ofav_site_ch ~ type, data = urban_ofav_site_ch)
summary(hazard_ofav_site_ch)

hazard_ofav_site_lh <- coxph(surv_ofav_site_lh ~ type, data = urban_ofav_site_lh)
summary(hazard_ofav_site_lh)

hazard_ssid_site_ch <- coxph(surv_ssid_site_ch ~ type, data = urban_ssid_site_ch)
summary(hazard_ssid_site_ch)

hazard_ssid_site_lh <- coxph(surv_ssid_site_lh ~ type, data = urban_ssid_site_lh)
summary(hazard_ssid_site_lh)

# hazard ratio plots
hazplot_ofav_site_ch <- ggforest(hazard_ofav_site_ch, data = urban_ofav_site_ch)
hazplot_ofav_site_ch

hazplot_ofav_site_lh <- ggforest(hazard_ofav_site_lh, data = urban_ofav_site_lh) + coord_flip(ylim=c(0,1))
hazplot_ofav_site_lh

hazplot_ssid_site_ch <- ggforest(hazard_ssid_site_ch, data = urban_ssid_site_ch)
hazplot_ssid_site_ch

hazplot_ssid_site_lh <- ggforest(hazard_ssid_site_lh, data = urban_ssid_site_lh)
hazplot_ssid_site_lh

# survivorship by treatment/type figure
urban_ofav_type_multiplot<-ggarrange(survival_ofav_site_ch$plot,
                                     survival_ofav_site_lh$plot, 
                                     survival_ofav_site_ch$table, 
                                     survival_ofav_site_lh$table,
                                     hazplot_ofav_site_ch,
                                     hazplot_ofav_site_lh,
                                     heights = c(2, 0.5, 0.5),
                                     ncol = 2, nrow = 3)
urban_ofav_type_multiplot

ggsave("urban survivorship ofav site type.pdf", urban_ofav_type_multiplot, width=10, height=10,dpi = 300)

urban_ssid_type_multiplot<-ggarrange(survival_ssid_site_ch$plot,
                                     survival_ssid_site_lh$plot, 
                                     survival_ssid_site_ch$table, 
                                     survival_ssid_site_lh$table,
                                     hazplot_ssid_site_ch,
                                     hazplot_ssid_site_lh,
                                     heights = c(2, 0.5, 0.5),
                                     ncol = 2, nrow = 3)
urban_ssid_type_multiplot

ggsave("urban survivorship ssid site type.pdf", urban_ssid_type_multiplot, width=10, height=10,dpi = 300)
