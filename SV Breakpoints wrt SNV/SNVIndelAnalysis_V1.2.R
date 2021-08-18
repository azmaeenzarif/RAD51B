library(reshape2)
library(data.table)
library(dplyr)
library(plyr)
library(rvest)
library(tidyverse)
library(ggplot2)
library(fields)
library(lattice)
library(diagram)
library(tidyr)

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Question 1/')

#Loads all RAD51B mutations (in 16 files) from PCAWG dataset
mutation_list <- list.files("RAD51B_Mutations/", pattern = '\\.tsv$', full.names = TRUE)
all_mutations <- lapply(mutation_list, function(x){
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})

#Combines all files into one data frame
all_mutations <- bind_rows(all_mutations)

#Data frame of SV calls
SV_calls <- data.frame(read.csv('merged_1.6.1.csv', header = TRUE))

#Selects SV calls for only Chr 14 (where RAD51B normally located)
chr14 <- SV_calls[SV_calls$seqnames=='14',]

#To separate locus of each SNV/Indel
all_mutations$SNV_Indel_Locus <- all_mutations$Genomic.DNA.Change

all_mutations$SNV_Indel_Locus <- substring(all_mutations$SNV_Indel_Locus, 9)

all_mutations$SNV_Indel_Locus <- gsub("[^0-9]", "",all_mutations$SNV_Indel_Locus)
  
all_mutations$SNV_Indel_Locus <- data.frame(all_mutations$SNV_Indel_Locus,1)

#To separate each type of SNV/Indel
all_mutations$SNV_Indel_Mutation <- all_mutations$Genomic.DNA.Change

all_mutations$SNV_Indel_Mutation <- substring(all_mutations$SNV_Indel_Mutation, 17)




#Ranks SNV/Indel in terms of frequency; rank denotes colour in plot
all_mutations$SNV_Indel_Mutation_Frequency <- table(all_mutations$SNV_Indel_Mutation)[all_mutations$SNV_Indel_Mutation]


all_mutations$SNV_Indel_Mutation_Rank <- rank(all_mutations$SNV_Indel_Mutation_Frequency, ties.method='min')



#Plot of all pairs of SV breakpoints in colour-coordinated manner and shows where each type of SNVs/Indels occur in relation to them
plot(all_mutations$SNV_Indel_Locus, type = 'b', xlab = 'Chr14 Locus', ylab = '', col = all_mutations$SNV_Indel_Mutation_Rank, ylim = c(0.5,1.5), yaxt = 'n')

View(all_mutations)

set.seed(1)
for (i in 1:nrow(chr14)) {
  xline(chr14$start[i])
  xline(chr14$altpos[i])
}

pdf('RAD51B SNVs-Indels wrt Chr 14 breakpoints.pdf', width = 30, height = 10)
dev.off()



