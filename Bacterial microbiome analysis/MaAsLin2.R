setwd ("/path/to/samples/")


library(Maaslin2)
library(edgeR)
library(metagenomeSeq)
library(pscl)
library(pbapply)
library(dplyr)
library(vegan)
library(chemometrics)
library(pheatmap)
library(cplm)
library(logging)
library(data.table)
library(lmerTest)

## Transpose relative abundance so samples are rows #
data<-read.delim("krakenUniq_clean_species_default_relative_abundance.txt", row.names = 1)
data.transposed<-t(data)
write.table(data.transposed, file="species_relative_abundance_trans.txt", sep="\t", row.names=TRUE)

###All
df_input_data = read.table("species_relative_abundance_trans.txt", header = TRUE, sep = "\t",
                           row.names = 1,
                           stringsAsFactors = FALSE)

df_input_metadata = read.table("MasterMetadataFile.txt", header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
fit_data = Maaslin2(
  input_data = df_input_data, 
  input_metadata = df_input_metadata, 
  normalization = "NONE",
  output = "results_personcode", 
  fixed_effects = c("Person_Code", "Timepoint", "InfectionStatusOG", "EverInfection", "HIVstatus", "Abx10days" ),
  reference = c("InfectionStatusOG","Negative","EverInfection","No","Timepoint", "T1"), #if you have more than two values in a column 
  random_effects = c("Study_ID"),
  standardize = TRUE)

#plot PersonCode
library(ggplot2)
library(RColorBrewer)
data<-read.delim("bubbleplot_personcode.txt")
ggplot(data, aes(x = PersonCode, y = SpeciesName, size = Prevalence, color = Mean_abundance)) +
  geom_point()+
  scale_size(range = c(0,7),
             breaks = c(0, 0.25, 0.5, 0.75, 1),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D"),
                        values = c(0, 0.05, 0.125, 0.5, 1))+
  theme_classic() +
  theme(axis.text.y = element_text(size = 3))  # Adjust font size as needed

# Mom/Infant SHAP
library(ggplot2)
library(RColorBrewer)
data<-read.delim("mom_baby_shap_count_mean.txt")
ggplot(data, aes(x = Dataset, y = OTU, size = Prevalence, color = Mean_ABS_SHAP_Value)) +
  geom_point()+
  scale_size(range = c(0,8),
             breaks = c(0, 0.25, 0.5, 0.75, 1),
             labels = c("0", "25", "50", "75", "100"))+
  scale_color_gradientn(colors = c("#9ECAE1", "#2171B5", "#EF3B2C", "#A50F15", "#67000D")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8))  # Adjust font size as needed