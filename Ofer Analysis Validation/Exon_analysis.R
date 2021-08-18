library(reshape)
library(data.table)
library(dplyr)
library(plyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(R.utils)
library(gridExtra)
library(tibble)
library(httr)
library(jsonlite)

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Validation Analysis/')


exon_coordinates <- data.frame(read.csv('RAD51B_Exon_Coordinates.csv'))


filtered_list <- list.files("RAD51B_Exon_Mutations/", pattern = '\\.tsv$', full.names = TRUE)
all_filtered_mutations <- lapply(filtered_list, function(x){
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})

exon_mutations <- bind_rows(all_filtered_mutations)

exon_mutations_single_base <- subset(exon_mutations, exon_mutations$Type == 'single base substitution')
#View(exon_mutations_single_base)
exon_mutations_insertion <- subset(exon_mutations, exon_mutations$Type == 'insertion of <=200bp')
#View(exon_mutations_insertion)
exon_mutations_del <- subset(exon_mutations, exon_mutations$Type == 'deletion of <=200bp')
#View(exon_mutations_del)
exon_mutations_multiple_base <- subset(exon_mutations, exon_mutations$Type == 'multiple base substitution (>=2bp and <=200bp)')
#View(exon_mutations_multiple_base)


exon_mutations_single_base$SNV_Indel_Locus <- exon_mutations_single_base$Genomic.DNA.Change

exon_mutations_single_base$SNV_Indel_Locus <- substring(exon_mutations_single_base$SNV_Indel_Locus, 9)

exon_mutations_single_base$SNV_Indel_Locus <- gsub("[^0-9]", "",exon_mutations_single_base$SNV_Indel_Locus)


exon_coordinates$Start <- gsub(',','',exon_coordinates$Start)
exon_coordinates$End <- gsub(',','',exon_coordinates$End)


transform(exon_coordinates, Start = as.numeric(Start), End = as.numeric(End))

exon_mutations_single_base$SNV_Indel_Locus <- gsub(',','',exon_mutations_single_base$SNV_Indel_Locus)

transform(exon_mutations_single_base, SNV_Indel_Locus = as.numeric(SNV_Indel_Locus))



exon_mutations_single_base$exon_status <- ifelse(sapply(exon_mutations_single_base$SNV_Indel_Locus, function(x) 
                                                any(x >= exon_coordinates$Start & x <= exon_coordinates$End)),'Exon','Intron')


#exon_mutations_single_base_final <- exon_mutations_single_base[ (exon_mutations_single_base$exon_status == 'Exon') | (exon_mutations_single_base[grep('Exon', exon_mutations_single_base$Consequences),])]
#subset(exon_mutations_single_base, (exon_mutations_single_base$exon_status == 'Exon') | (exon_mutations_single_base[grep('Exon', exon_mutations_single_base$Consequences),]) )


exon_mutations_single_base_final <- bind_rows((subset(exon_mutations_single_base, (exon_mutations_single_base$exon_status == 'Exon'))), (subset(exon_mutations_single_base[grep('Exon', exon_mutations_single_base$Consequences),])))
exon_mutations_single_base_final <- exon_mutations_single_base_final[!duplicated(exon_mutations_single_base_final),]




exon_mutations_del_final <- subset(exon_mutations_del[grep('Exon', exon_mutations_del$Consequences),])
exon_mutations_insertion_final <- subset(exon_mutations_insertion[grep('Exon', exon_mutations_insertion$Consequences),])
exon_mutations_multiple_base_final <- subset(exon_mutations_multiple_base[grep('Exon', exon_mutations_multiple_base$Consequences),])


exon_mutations_final <- bind_rows(exon_mutations_single_base_final, exon_mutations_multiple_base_final, exon_mutations_insertion_final, exon_mutations_del_final)

View(exon_mutations_final)
write.csv(exon_mutations_final, 'RAD51B_exon_mutations_final.csv')

