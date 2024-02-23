setwd("/path/to/samples")

library(vegan)
library(ggplot2)
library(ape)
library(dplyr)
library(reshape2)
library(phyloseq)
library(nlme)
library(RColorBrewer)

#Read in files 
set.seed(100)
feature_table<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("MasterMetadataFile.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
write.table(axes, "BC_PCoA_axes_all_samples.txt", sep = "\t")

data<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
write.table(dis.bc2,"BC.txt", sep = '\t') #Matrix with double values
dis.bc2[lower.tri(dis.bc2, diag = TRUE)] <- NA #Removes double values
dis.bc2.long<-melt(dis.bc2)
write.table(dis.bc2.long, "BC_long_format_102123.txt", sep="\t", row.names = FALSE)

#Mother/Infant PCoA and PERMANOVA
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "Person_Code", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("Infant","Mother"),
                     values = c("darkgreen","deeppink")) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

identical(metadata$I, rownames(dis.bc2))
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"))
adonis.mom.infant <- adonis2(dis.bc ~ Person_Code + Timepoint + EverInfected + HIVstatus + InfectionStatusOG + Covid_Status + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.mom.infant

#MotherHyp1
metadata<-import_qiime_sample_data("MasterMetaMom.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
set.seed(100)
physeq_ctrls <- subset_samples(physeq, Person_Code=="Mother")
ord_ctrls <- ordinate(physeq_ctrls,
                      method = "PCoA",
                      distance = "bray")
axes.ctrls <- as.matrix(ord_ctrls$vectors)
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "Covid_Status", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("NeverPositive","Negative","Positive"),
                     values = c("gray","lightpink","deeppink")) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

data<-read.delim("krakenUniq_clean_species_default_relative_abundance_mom.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("MasterMetaMom.txt", header = TRUE, sep = '\t')
identical(metadata$I, rownames(dis.bc2))

set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Study_ID)) ##Control for repeated sampling. Replace Patient_factor with appropriate variable from your metadata.
adonis.Hyp1M <- adonis2(dis.bc ~ Covid_Status + Timepoint + EverCovid +HIVstatus + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.Hyp1M

#InfantHyp1
set.seed(100)
feature_table<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("MasterMetaInfant.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
set.seed(100)
physeq_ctrls <- subset_samples(physeq, Person_Code=="Infant")
ord_ctrls <- ordinate(physeq_ctrls,
                      method = "PCoA",
                      distance = "bray")
axes.ctrls <- as.matrix(ord_ctrls$vectors)
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "Infection_Status", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("NeverPositive","Negative","Positive","RVP"),
                     values = c("gray","lightgreen","darkgreen","purple")) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

data<-read.delim("krakenUniq_clean_species_default_relative_abundance_baby.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("MasterMetaInfant.txt", header = TRUE, sep = '\t')
identical(metadata$I, rownames(dis.bc2))

set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Study_ID)) ##Control for repeated sampling. Replace Patient_factor with appropriate variable from your metadata.
adonis.Hyp1B <- adonis2(dis.bc ~ Infection_Status + Timepoint + EverCovid + HIVstatus + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.Hyp1B

#MotherHyp2
set.seed(100)
feature_table<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("LME_Mothers_FINAL.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
set.seed(100)
physeq_ctrls <- subset_samples(physeq, Person_Code=="Mother")
ord_ctrls <- ordinate(physeq_ctrls,
                      method = "PCoA",
                      distance = "bray")
axes.ctrls <- as.matrix(ord_ctrls$vectors)
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "Covid_Status2", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("Negative","Positive","Recovery"),
                     values = c("gray","lightpink","deeppink")) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

data<-read.delim("krakenUniq_clean_species_default_relative_abundance_Hyp2M.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("LME_Mothers_FINAL.txt", header = TRUE, sep = '\t')
identical(metadata$I, rownames(dis.bc2))
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Study_ID)) ##Control for repeated sampling. Replace Patient_factor with appropriate variable from your metadata.
adonis.Hyp2B <- adonis2(dis.bc ~ Covid_Status2 + Timepoint + EverCovid + HIVStatus + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.Hyp2B


#Mothers HIV
set.seed(100)
feature_table<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("Mothers_FINAL.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
set.seed(100)
physeq_ctrls <- subset_samples(physeq, Person_Code=="Mother")
ord_ctrls <- ordinate(physeq_ctrls,
                      method = "PCoA",
                      distance = "bray")
axes.ctrls <- as.matrix(ord_ctrls$vectors)
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "HIVStatus", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c('#4575b4','#d53e4f')) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

data<-read.delim("krakenUniq_clean_species_default_relative_abundance_Hyp2M.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("LME_Mothers_FINAL.txt", header = TRUE, sep = '\t')
identical(metadata$I, rownames(dis.bc2))
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Study_ID)) ##Control for repeated sampling. Replace Patient_factor with appropriate variable from your metadata.
adonis.MomHIV <- adonis2(dis.bc ~ Covid_Status2 + Timepoint + EverCovid + HIVStatus + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.MomHIV


#InfantHyp2
set.seed(100)
feature_table<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("Infant_FINAL.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
set.seed(100)
physeq_ctrls <- subset_samples(physeq, Person_Code=="Infant")
ord_ctrls <- ordinate(physeq_ctrls,
                      method = "PCoA",
                      distance = "bray")
axes.ctrls <- as.matrix(ord_ctrls$vectors)
#write.table(axes.ctrls, "BC_PCoA_axes_ctrls.txt", sep = "\t")
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "HIVStatus", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("Negative","Positive"),
                     values = c('#4575b4','#d53e4f')) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

data<-read.delim("krakenUniq_clean_species_default_relative_abundance_Hyp2I.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("LME_Infant_FINAL.txt", header = TRUE, sep = '\t')
identical(metadata$I, rownames(dis.bc2))
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Study_ID)) ##Control for repeated sampling. Replace Patient_factor with appropriate variable from your metadata.
adonis.BabyHIV <- adonis2(dis.bc ~ Infection_Status2 + Timepoint + EverCovid + HIVStatus + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.BabyHIV

#InfantHyp2
set.seed(100)
feature_table<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
feature_table <- as.data.frame(feature_table)
dim(feature_table)
OTU=otu_table(feature_table,taxa_are_rows =TRUE)
sample_names(OTU)
metadata<-import_qiime_sample_data("Infant_FINAL.txt")
sample_names(metadata)
physeq =             phyloseq(OTU,metadata)
ord <- ordinate(physeq,
                method = "PCoA",
                distance = "bray")
axes <- as.matrix(ord$vectors)
set.seed(100)
physeq_ctrls <- subset_samples(physeq, Person_Code=="Infant")
ord_ctrls <- ordinate(physeq_ctrls,
                      method = "PCoA",
                      distance = "bray")
axes.ctrls <- as.matrix(ord_ctrls$vectors)
PCoA.BC  <- plot_ordination(physeq = physeq,
                            ordination = ord,
                            #shape = "Person_Code", # metadata variable
                            color = "HIVStatus", # metadata variable
                            axes = c(1,2),
                            title='Bray Curtis PCoA') +
  stat_ellipse() +
  theme_bw() +
  scale_color_manual(breaks = c("Positive","Negative"),
                     values = c("gray","lightgreen","darkgreen","purple","darkblue")) +
  theme(text=element_text(size=20))+
  geom_point(size = 6)
PCoA.BC

data<-read.delim("krakenUniq_clean_species_default_relative_abundance_Hyp2I.txt", row.names = 1)
data.transposed<-t(data)
dis.bc <- vegdist(data.transposed, method = "bray")
dis.bc2<-as.matrix(dis.bc)
metadata <- read.delim("LME_Infant_FINAL.txt", header = TRUE, sep = '\t')
identical(metadata$I, rownames(dis.bc2))
set.seed(100)
perm <- how(nperm = 999, within = Within(type = "free"), plots = Plots(strata = metadata$Study_ID)) ##Control for repeated sampling. Replace Patient_factor with appropriate variable from your metadata.
adonis.Hyp2B <- adonis2(dis.bc ~ Infection_Status2 + Timepoint + EverCovid + HIVStatus + Abx10days, data = metadata, permutations = perm, by = "margin")
adonis.Hyp2B
