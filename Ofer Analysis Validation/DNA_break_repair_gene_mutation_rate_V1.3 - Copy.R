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




RAD51B_Mutated_Exon_SV <- data.frame(fread('RAD51B_Mutated_Exon_SV.gz'))

RAD51B_Mutated_Exon_CNV <- data.frame(fread('RAD51B_Mutated_Exon_CNV.gz'))

View(RAD51B_Mutated_Exon_SV)
View(RAD51B_Mutated_Exon_CNV)

length(unique(RAD51B_Mutated_Exon_SV$icgc_donor_id)) #190
length(unique(RAD51B_Mutated_Exon_SV$icgc_sample_id)) #240
length(unique(RAD51B_Mutated_Exon_SV$sv_id)) #35829




RAD51B_Mutated_Exon_SV$icgc_sample_id <- substring(RAD51B_Mutated_Exon_SV$icgc_sample_id, 3)

#RAD51B_Not_Mutated_SV[order(RAD51B_Not_Mutated_SV$icgc_sample_id),]
RAD51B_Mutated_Exon_SV[order(RAD51B_Mutated_Exon_SV$icgc_sample_id),]



#RAD51B_Not_Mutated_SV_per_Sample <- data.table(RAD51B_Not_Mutated_SV)[,.N, by=RAD51B_Not_Mutated_SV$icgc_sample_id]
#setnames(RAD51B_Not_Mutated_SV_per_Sample, 'N', 'RAD51B WT SV Mutation Frequency')


RAD51B_Mutated_Exon_SV_per_Sample <- data.table(RAD51B_Mutated_Exon_SV)[,.N, by=RAD51B_Mutated_Exon_SV$icgc_sample_id]
setnames(RAD51B_Mutated_Exon_SV_per_Sample, 'N', 'RAD51B Mutant SV Mutation Frequency')


View(RAD51B_Mutated_Exon_SV_per_Sample)

mean(RAD51B_Mutated_Exon_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`) #162.8125
median(RAD51B_Mutated_Exon_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`) #80






#RAD51B_gene_SV <- RAD51B_Mutated_Exon_SV[(RAD51B_Mutated_Exon_SV$chr_from =='14' & RAD51B_Mutated_Exon_SV$chr_from_bkpt>='68286496' & RAD51B_Mutated_Exon_SV$chr_from_bkpt<= '69196935') | (RAD51B_Mutated_Exon_SV$chr_to =='14' & RAD51B_Mutated_Exon_SV$chr_to_bkpt>='68286496' & RAD51B_Mutated_Exon_SV$chr_to_bkpt<= '69196935'),]

#View(RAD51B_gene_SV)






#Count of number of donors with germline data for analysis

RAD51B_Mutated_Exon_SV_donors <- distinct(as.data.frame(RAD51B_Mutated_Exon_SV$icgc_donor_id))


SGV_count <- 0


for (i in 1:nrow(RAD51B_Mutated_Exon_SV_donors)) {
  donors <- toString(RAD51B_Mutated_Exon_SV_donors[i,])
  get_donors <- GET(gsub(' ', '', paste('https://dcc.icgc.org/api/v1/donors/',donors)))
  print(get_donors)
  get_donors_text <- content(get_donors, 'text', encoding ='UTF-8')
  print(get_donors_text)
  get_donors_json <- fromJSON(get_donors_text, flatten = TRUE)
  print(get_donors_json)
  get_donors_json_df <- data.frame(get_donors_json = unlist(get_donors_json))
  print(get_donors_json)

 
  if (get_donors_json_df['sgvExists',] == TRUE) {
    SGV_count <- SGV_count+1
    
  }
}


SGV_count #= 36 donors with germline data to analyse (1/5th)


#Exon Analysis

exon_list <- list.files("RAD51B_Exon_Mutations/", pattern = '\\.tsv$', full.names = TRUE)
all_exon_mutations <- lapply(exon_list, function(x){
  read.table(file = x,
             sep = '\t',
             header = TRUE)
  
})

exon_mutations_initial <- bind_rows(all_exon_mutations)


distinct(as.data.frame(exon_mutations_initial$Type))


exon_mutations_final <- exon_mutations_initial[-grep('Intron', exon_mutations_initial$Consequences),]


#write.csv(exon_mutations_final,'RAD51B_exon_mutations_final.csv')







