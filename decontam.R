##Contamination analysis
setwd('/path/to/samples/')

##Decontam on krakenUniq bacterial species
#install.packages("BiocManager")
#BiocManager::install("decontam")
library(decontam)

##Method: prevalence. Include (+) control. Use this.
data <- read.delim("masked_combined_species.txt", row.names=1)
metadata<-read.delim("Metadata_decontam.txt", row.names = 1)
data<-as.matrix(t(data))

contam.default<-isContaminant(data, method = 'prevalence', neg =metadata$neg_ctrl, threshold=0.1)
write.csv(contam.default, "krakenUniq_decontam_default.csv")
contam.intermediate<-isContaminant(data, method = 'prevalence', neg =metadata$neg_ctrl, threshold=0.25)
write.csv(contam.intermediate, "krakenUniq_decontam_intermediate.csv")
contam.strict<-isContaminant(data, method = 'prevalence', neg =metadata$neg_ctrl, threshold=0.5)
write.csv(contam.strict, "krakenUniq_decontam_strict.csv")

clean.default <- t(data[,!contam.default$contaminant=="TRUE"])
clean.default <- clean.default[,!metadata$neg_ctrl=="TRUE"]
clean.default <- clean.default[rowSums(clean.default) > 0,]
write.table(clean.default,"krakenUniq_clean_species_default.txt", sep="\t", row.names = TRUE, col.names = NA)

clean.intermediate <- t(data[,!contam.intermediate$contaminant=="TRUE"])
clean.intermediate <- clean.intermediate[,!metadata$any_ctrl=="TRUE"]
clean.intermediate <- clean.intermediate[rowSums(clean.intermediate) > 0,]
write.table(clean.intermediate,"krakenUniq_clean_species_intermediate.txt", sep="\t", row.names = TRUE, col.names = NA)

clean.strict <- t(data[,!contam.strict$contaminant=="TRUE"])
clean.strict <- clean.strict[,!metadata$any_ctrl=="TRUE"]
clean.strict <- clean.strict[rowSums(clean.strict) > 0,]
write.table(clean.strict,"krakenUniq_clean_species_strict.txt", sep="\t", row.names = TRUE, col.names = NA)