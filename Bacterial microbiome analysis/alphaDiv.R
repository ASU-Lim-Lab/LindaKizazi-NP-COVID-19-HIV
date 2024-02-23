setwd("/path/to/samples")

library(vegan)
library(ggplot2)
library(nlme)
library (RColorBrewer)

## General 
data<-read.delim("NP_subset_mother_infant", row.names = 1)
dataTransposed<-t(data)
shannon<-diversity(dataTransposed, index = 'shannon')
write.table(shannon,"NP_shannon_012323.txt", sep = '\t')

##Mother/Infant plots and LME
plot.data <- read.delim("alphadiv.txt")
richness.plot <- ggplot(plot.data, aes(x = Timepoint, y = Richness, color = Person_Code, fill = Person_Code)) +
  stat_smooth(method="loess", se=TRUE, span=0.75, level=0.95) +
  geom_point() +
  scale_x_continuous(breaks=c(1,2,3,4)) +
  theme_bw()
richness.plot

plot.data <- read.delim("alphadiv.txt")
shannon.plot <- ggplot(plot.data, aes(x = Timepoint, y = Shannon, color = Person_Code, fill = Person_Code)) +
  stat_smooth(method="loess", se=TRUE, span=0.75, level=0.95) +
  geom_point() +
  scale_x_continuous(breaks=c(1,2,3,4)) +
  theme_bw()
shannon.plot

richness.no.interaction <- lme(Richness ~ Timepoint + EverCovid + CovidStatusOG + Person_Code + HIVStatus + Abx10days, random =~1|Patient, data = plot.data)
summary(richness.no.interaction)

shannon.no.interaction <- lme(Shannon ~ Timepoint + EverCovid + CovidStatusOG + Person_Code + HIVStatus + Abx10days, random =~1|Patient, data = plot.data)
summary(shannon.no.interaction)

#Richness over time no interaction term for single relationship.#
richness.no.interaction <- lme(Richness ~ Person_Code, random =~1|Patient, data = plot.data)
richness_summary <- summary(richness.no.interaction)
richness_summary$tTable

#Shannon index over time no interaction term for single relationship.#
shannon.no.interaction <- lme(Shannon ~ Person_Code, random =~1|Patient, data = plot.data)
shannon_summary <- summary(shannon.no.interaction)
shannon_summary$tTable

pvalues.5 <- c(0.1030, 0.7756) #alpha
p.adjust(pvalues.5,method="BH")
pvalues.6 <- c(0.5368, 0.7789) #alpha
p.adjust(pvalues.6,method = "BH")
pvalues.qpcr5 <- c(0.3755, 0.2408)
p.adjust(pvalues.qpcr5,method="BH")

## Mothers HIV
plot.data <- read.delim("LME_Mothers.txt")
shannon.plot <- ggplot(plot.data, aes(x = Timepoint, y = Shannon, color = HIVStatus, fill = HIVStatus)) +
  stat_smooth(method="loess", se=TRUE, span=0.75, level=0.95) +
  scale_shape_manual(values=c("Positive", "Negative")) + 
  scale_color_manual(values=c('#4575b4','#d53e4f'))+
  scale_fill_manual(values=c('#4575b4','#d53e4f')) +
  geom_point() +
  scale_x_continuous(breaks=c(1,2,3,4)) +
  theme_bw()
shannon.plot

richness.no.interaction <- lme(Richness ~ Timepoint + EverCovid + CovidStatusOG + Abx10days + HIVStatus, random =~1|Patient, data = plot.data)
summary(richness.no.interaction)

shannon.no.interaction <- lme(Shannon ~ Timepoint + EverCovid + CovidStatusOG + Abx10days + HIVStatus, random =~1|Patient, data = plot.data)
summary(shannon.no.interaction)

## Infants HIV
plot.data <- read.delim("LME_Infant.txt")
shannon.plot <- ggplot(plot.data, aes(x = Timepoint, y = Shannon, color = HIVStatus, fill = HIVStatus)) +
  stat_smooth(method="loess", se=TRUE, span=0.75, level=0.95) +
  scale_shape_manual(values=c("Positive", "Negative")) + 
  scale_color_manual(values=c('#4575b4','#d53e4f'))+
  scale_fill_manual(values=c('#4575b4','#d53e4f')) +
  geom_point() +
  scale_x_continuous(breaks=c(1,2,3,4)) +
  theme_bw()
shannon.plot

richness.no.interaction <- lme(Richness ~ Timepoint + EverInfection + InfectionStatusOG + Abx10days + HIVStatus, random =~1|Patient, data = plot.data)
summary(richness.no.interaction)

shannon.no.interaction <- lme(Shannon ~ Timepoint + EverInfection + InfectionStatusOG + Abx10days + HIVStatus, random =~1|Patient, data = plot.data)
summary(shannon.no.interaction)


