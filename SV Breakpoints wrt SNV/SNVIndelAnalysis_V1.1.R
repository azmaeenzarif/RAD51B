library('httr')
library('jsonlite')
library(reshape2)
library(data.table)
library(dplyr)
library(plyr)
library(rvest)
library(tidyverse)
library(ggplot2)
library(fields)
library(lattice)

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Question 1/')

# MutID <- as.data.frame(fread('mutation-ids-for-set-RAD51B, Chr14_67819779-68730218, PCAWG.tsv', header = FALSE))
# 
# View(MutID)
# nrow(MutID)
# 
# get_mutations <- GET('https://dcc.icgc.org/api/v1/mutations/MU30396090')
# 
# get_mutations
# 
# get_mutations_text <- content(get_mutations, 'text')
# 
# get_mutations_text
# 
# get_mutations_json <- fromJSON(get_mutations_text, flatten = TRUE)
# 
# 
# 
# get_mutations_json_df <- data.frame(get_mutations_json = unlist(get_mutations_json))
# 
# nrow(get_mutations_json_df)
# 
# mutations_df$first <- get_mutations_json_df
# 
# View(mutations_df)

# 
# 
# get_mutations1 <- GET('https://dcc.icgc.org/api/v1/mutations/MU32937387')
# 
# 
# get_mutations_text1 <- content(get_mutations1, 'text')
# 
# 
# get_mutations_json1 <- fromJSON(get_mutations_text1, flatten = TRUE)
# 
# get_mutations_json1_df <- data.frame(get_mutations_json1 = unlist(get_mutations_json1))
# 
# mutations_df$alt <- get_mutations_json1_df




#Create empty dataframe to store mutations
# mutations_df <- data.frame(matrix(nrow = 87))
# 
# View(mutations_df)


#Function to parse mutations using mutation IDs of PCAWG RAD51B SNVs/Indels
# SNVIndelSummary <- function(mutation_id, df, i) {
#   View(df)
#   mut <- toString(mutation_id)
#   print (mut)
#   get_mutation <- GET(gsub(' ', '', paste('https://dcc.icgc.org/api/v1/mutations/',mut)))
#   get_mutation_text <- content(get_mutation, 'text', encoding ='UTF-8')
#   print(get_mutation_text)
#   get_mutation_json <- fromJSON(get_mutation_text, flatten = TRUE)
#   print(get_mutation_json)
#   get_mutation_json_df <- data.frame(get_mutation_json = unlist(get_mutation_json))
#   print(get_mutation_json_df)
#   df <- merge(df, get_mutation_json_df, all=TRUE)
#   View(df)
# }


#For-loop to load every mutation and extract relevant data into mutations_df dataframe
# for i in xrange(0,nrow(MutID),batchsize) {
#   row <- MutID[i,]
#   print(row)
#   SNVIndelSummary(row, mutations_df, i)
# }



mutation_list <- list.files("RAD51B_Mutations/", pattern = '\\.tsv$', full.names = TRUE)
all_mutations <- lapply(mutation_list, function(x){
  read.table(file = x,
             sep = '\t',
             header = TRUE)
  
  
})


all_mutations <- bind_rows(all_mutations)

View(all_mutations)


SV_calls <- data.frame(read.csv('merged_1.6.1.csv', header = TRUE))

View(SV_calls)

colnames(SV_calls)

chr14 <- SV_calls[SV_calls$seqnames=='14',]


all_mutations$SNV_Indel_Locus <- all_mutations$Genomic.DNA.Change

all_mutations$SNV_Indel_Locus <- substring(all_mutations$SNV_Indel_Locus, 9)

all_mutations$SNV_Indel_Locus <- gsub("[^0-9]", "",all_mutations$SNV_Indel_Locus)
  
all_mutations$SNV_Indel_Locus <- data.frame(all_mutations$SNV_Indel_Locus,1)


 

all_mutations$SNV_Indel_Mutation <- all_mutations$Genomic.DNA.Change



all_mutations$SNV_Indel_Mutation <- substring(all_mutations$SNV_Indel_Mutation, 17)

all_mutations$SNV_Indel_Mutation_Frequency <- table(all_mutations$SNV_Indel_Mutation)[all_mutations$SNV_Indel_Mutation]


all_mutations$SNV_Indel_Mutation_Rank <- rank(all_mutations$SNV_Indel_Mutation_Frequency, ties.method='min')



View(all_mutations)



plot(all_mutations$SNV_Indel_Locus, type = 'b', xlab = 'Chr14 Locus', ylab = '', yaxt = 'n', col = all_mutations$SNV_Indel_Mutation_Rank)

for (i in 1:nrow(chr14)) {
  xline(chr14$start[i], col = i)
  xline(chr14$altpos[i], col = i)
}





