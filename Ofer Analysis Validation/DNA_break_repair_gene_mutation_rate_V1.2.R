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



#Using merged_1.6.1.csv file


SV_calls <- data.frame(read.csv('merged_1_6_1_w_ofer_annot.csv', header = TRUE))


#To load RAD51B SV loci, once set, csv file with donor ids (manually entered) now used 

#chr14 <- SV_calls[SV_calls$seqnames=='14',]

#RAD51B_locus <- chr14[chr14$start>='68286496' & chr14$start<='69196935' & chr14$altpos>='68286496' & chr14$altpos<='69196935' ,]


RAD51B_locus <- data.frame(read.csv('RAD51B_locus.csv', header = TRUE))

View(RAD51B_locus)

RAD51B_Not_Mutated_SV <- data.frame(fread('RAD51B_Not_Mutated_SV.gz'))

RAD51B_Mutated_SV <- data.frame(fread('RAD51B_Mutated_SV.gz'))

RAD51B_Not_Mutated_SV$icgc_sample_id <- substring(RAD51B_Not_Mutated_SV$icgc_sample_id, 3)
RAD51B_Mutated_SV$icgc_sample_id <- substring(RAD51B_Mutated_SV$icgc_sample_id, 3)

RAD51B_Not_Mutated_SV[order(RAD51B_Not_Mutated_SV$icgc_sample_id),]
RAD51B_Mutated_SV[order(RAD51B_Mutated_SV$icgc_sample_id),]



RAD51B_Not_Mutated_SV_per_Sample <- data.table(RAD51B_Not_Mutated_SV)[,.N, by=RAD51B_Not_Mutated_SV$icgc_sample_id]
setnames(RAD51B_Not_Mutated_SV_per_Sample, 'N', 'RAD51B WT SV Mutation Frequency')


RAD51B_Mutated_SV_per_Sample <- data.table(RAD51B_Mutated_SV)[,.N, by=RAD51B_Mutated_SV$icgc_sample_id]
setnames(RAD51B_Mutated_SV_per_Sample, 'N', 'RAD51B Mutant SV Mutation Frequency')



mocklist <- list()

RAD51B_DO496 <- GET('https://dcc.icgc.org/api/v1/donors/DO496')
RAD51B_DO496_text <- content(RAD51B_DO496, 'text', encoding = 'UTF-8')
RAD51B_DO496_json <- fromJSON(RAD51B_DO496_text, flatten=TRUE)
RAD51B_DO496_json_df <- data.frame(RAD51B_DO496_json = unlist(RAD51B_DO496_json))

RAD51B_DO41398 <- GET('https://dcc.icgc.org/api/v1/donors/DO41398')
RAD51B_DO41398_text <- content(RAD51B_DO41398, 'text', encoding = 'UTF-8')
RAD51B_DO41398_json <- fromJSON(RAD51B_DO41398_text, flatten=TRUE)
RAD51B_DO41398_json_df <- data.frame(RAD51B_DO41398_json = unlist(RAD51B_DO41398_json))

mocklist[[1]] <- RAD51B_DO496_json_df
mocklist[[2]] <- RAD51B_DO41398_json_df

View(RAD51B_DO41398_json_df['sgvExists',])


(matrix(unlist(mocklist), nrow=40))




RAD51B_donor_summary <- list()

SGV_count <- 0
  


for (i in 1:nrow(RAD51B_locus)) {
  donors <- toString(RAD51B_locus$donor_id[i])
  get_donors <- GET(gsub(' ', '', paste('https://dcc.icgc.org/api/v1/donors/',donors)))
  get_donors_text <- content(get_donors, 'text', encoding ='UTF-8')
  get_donors_json <- fromJSON(get_donors_text, flatten = TRUE)
  get_donors_json_df <- data.frame(get_donors_json = unlist(get_donors_json))
  RAD51B_donor_summary[[i]] <- get_donors_json_df
  
  if (get_donors_json_df['sgvExists',] == TRUE) {
    SGV_count <- SGV_count+1
  }
}








  
