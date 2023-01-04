#### packages ####

library(tidyverse)
library(dplyr)
library(ggpubr)

#### temp data import ####

templog <- read.csv(file="urban_stress_temp.csv", head=T) 
templog$Time<- lubridate::mdy_hms(templog$Time) # reformats time as month day year hours minutes seconds
templog$OA <- as.factor(templog$OA)
templog$Thermal <- as.factor(templog$Thermal)
str(templog)

temp_melt <- gather(templog, spp.ph.temp, measurement, Ssid.C.C:Ofav2.L.H, factor_key=TRUE) # melts dataframe to wide format
temp_melt %>%
  separate(spp.ph.temp, into=c("species","ph","temp"), sep="[.]") -> temp_melt # separates species and ph/temp factors to columns

temp_melt$species <-  as.factor(temp_melt$species) 
temp_melt$ph <-  as.factor(temp_melt$ph) 
temp_melt$temp <-  as.factor(temp_melt$temp) 

temp_melt$hour <-lubridate::hour(temp_melt$Time) # creates column of simplified hour by time


#### temp plot ####

fill.color<-c("#92c5de","#f4a582") # custom palette

temp_melt %>% # plot of temp means by hour
  filter(Thermal =="thermal") %>%
  group_by(hour, temp)  %>% 
  summarise(avg_temp = mean(measurement, na.rm=TRUE), sd = sd(measurement, na.rm=TRUE), n = n(), se = sd/sqrt(n)) %>% 
  ggplot(., aes(x=hour, y=avg_temp, color = temp)) + 
  theme_classic() +
  geom_ribbon(aes(fill=temp, group = temp, ymin = avg_temp - sd, ymax = avg_temp + sd), alpha = 0.2) +
  geom_line(aes(group = temp)) +
  scale_color_manual(values = fill.color) +
  scale_fill_manual(values = fill.color) +
  labs(y="Temperature (C)") -> temp_plot
plot(temp_plot)

temp_melt %>% 
  filter(Thermal =="thermal") %>%
  group_by(temp)  %>% 
  summarise(avg_temp = mean(measurement, na.rm=TRUE), sd = sd(measurement, na.rm=TRUE), n = n(), se = sd/sqrt(n)) -> temp_mean # table of treatment means
write.csv(temp_mean, file = "temp_mean.csv")

#### ph data import ####

phlog <- read.csv(file="urban_stress_ph.csv", head=T) 
phlog$Time<- lubridate::mdy_hms(phlog$Time)
phlog$Timepoint <- as.factor(phlog$Timepoint)
str(phlog)

ph_melt <- gather(phlog, spp.ph.temp, measurement, Ssid.C.C:Ofav.L.H, factor_key=TRUE)
ph_melt %>%
  separate(spp.ph.temp, into=c("species","ph","temp"), sep="[.]") -> ph_melt

ph_melt$species <-  as.factor(ph_melt$species) 
ph_melt$ph <-  as.factor(ph_melt$ph) 
ph_melt$temp <-  as.factor(ph_melt$temp) 

ph_melt$hour <-lubridate::hour(ph_melt$Time)


#### ph plot ####

fill.color2<-c("#92c5de","#0571b0")

ph_melt %>% 
  filter(OA =="OA") %>%
  group_by(hour, ph) %>% 
  summarise(avg_ph = mean(measurement, na.rm=TRUE), sd = sd(measurement, na.rm=TRUE), n = n(), se = sd/sqrt(n)) %>% 
  ggplot(., aes(x=hour, y=avg_ph, color = ph)) + 
  theme_classic() +
  geom_ribbon(aes(fill=ph, group = ph, ymin = avg_ph - sd, ymax = avg_ph + sd), alpha = 0.2) +
  geom_line(aes(group = ph)) +
  scale_color_manual(values = fill.color2) +
  scale_fill_manual(values = fill.color2) +
  labs(y="pH") -> ph_plot
plot(ph_plot)

ph_melt %>% 
  filter(OA =="OA") %>%
  group_by(hour, ph) %>% 
  summarise(avg_ph = mean(measurement, na.rm=TRUE), sd = sd(measurement, na.rm=TRUE), n = n(), se = sd/sqrt(n)) -> ph_mean_hour # table of treatment means by hour
write.csv(ph_mean_hour, file = "ph_mean_hour.csv")

ph_melt %>% 
  filter(OA =="OA") %>%
  group_by(ph) %>% 
  summarise(avg_ph = mean(measurement, na.rm=TRUE), sd = sd(measurement, na.rm=TRUE), n = n(), se = sd/sqrt(n)) -> ph_mean
write.csv(ph_mean, file = "ph_mean.csv")


#### multi plot ####

temp_table <- ggtexttable(temp_mean)

ph_table <- ggtexttable(ph_mean)
ph_table_hour <- ggtexttable(ph_mean_hour)


multiplot <- ggarrange(temp_plot, ph_plot, temp_table, ph_table, heights = c(3,1))
multiplot
ggsave("urban tank treatments.pdf", plot = multiplot, dpi = 300, width = 12, height = 6)
