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
library(Biostrings)


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

all_mutations$SNV_Indel_Locus <- as.integer(gsub("[^0-9]", "",all_mutations$SNV_Indel_Locus))
  

all_mutations$SNV_Indel_Chr <- as.integer(substring(all_mutations$Genomic.DNA.Change, 4,5))

View(all_mutations$SNV_Indel_Chr)

  
View(all_mutations_RAD51B_SNV_Indel)  

#To separate each type of SNV/Indel
all_mutations$SNV_Indel_Mutation <- all_mutations$Genomic.DNA.Change

all_mutations$SNV_Indel_Mutation <- substring(all_mutations$SNV_Indel_Mutation, 17)


#Ranks SNV/Indel in terms of frequency; rank denotes colour in plot
# all_mutations$SNV_Indel_Mutation_Frequency <- table(all_mutations$SNV_Indel_Mutation)[all_mutations$SNV_Indel_Mutation]
# 
# 
# all_mutations$SNV_Indel_Mutation_Rank <- rank(all_mutations$SNV_Indel_Mutation_Frequency, ties.method='min')



#Plot of all pairs of SV breakpoints in colour-coordinated manner and shows where each type of SNVs/Indels occur in relation to them
# plot(all_mutations$SNV_Indel_Locus, type = 'b', xlab = 'Chr14 Locus', ylab = '', col = all_mutations$SNV_Indel_Mutation_Rank, ylim = c(0.5,1.5), yaxt = 'n')
# 
# View(all_mutations)
# 
# set.seed(1)
# for (i in 1:nrow(chr14)) {
#   xline(chr14$start[i])
#   xline(chr14$altpos[i])
# }
# 
# pdf('RAD51B SNVs-Indels wrt Chr 14 breakpoints.pdf', width = 30, height = 10)
# dev.off()



exon_coordinates <- data.frame(read.csv('RAD51B_Exon_Coordinates.csv'))


exon_coordinates$Start <- gsub(',','',exon_coordinates$Start)
exon_coordinates$End <- gsub(',','',exon_coordinates$End)


transform(exon_coordinates, Start = as.numeric(Start), End = as.numeric(End))


all_mutations$exon_intron_status <- ifelse(sapply(all_mutations$SNV_Indel_Locus, function(x) 
  any(x >= exon_coordinates$Start & x <= exon_coordinates$End)),'Exon','Intron')




RAD51B_exon <- readDNAStringSet('exon_Homo_sapiens_RAD51B_sequence.fa')
RAD51B_intron <- readDNAStringSet('intron_Homo_sapiens_RAD51B_sequence.fa')


sum(nrow(as.data.frame(unlist(RAD51B_exon))), nrow(as.data.frame(unlist(RAD51B_intron))))

#totalbp of RAD51B gene = 910439
#total SNV/Indels = 1570
#If randomly distributed, would expect 1 SNV/Indel per 580 bp

random_exon_count_summary <- list()

random_intron_count_summary <- list()


# exon_intron_classifier <- function(x) {
#   
#   if (any(x >= exon_coordinates$Start & x <= exon_coordinates$End)) {
#     random_exon_count <<- random_exon_count + 1
#   } else {
#     random_intron_count <<- random_intron_count + 1 
#   }   
#   
#   
# }

# expected_random_loci_hits$exon_intron_status <- for (i in 1:nrow(expected_random_loci_hits)) {
#   ifelse(sapply(expected_random_loci_hits[i,], 
# }


# example <- as.integer(expected_random_loci_hits[1])
# typeof(example)
# exon_intron_classifier(example)
# random_intron_count

n <- 10000

for (y in 1:n) {
  
  expected_random_loci_hits <- sample(67819779:68730218, 1570, replace = TRUE)
  
  random_exon_count <- 0
  random_intron_count <- 0
  
  for (i in 1:length(expected_random_loci_hits)) {
    locus <- as.integer(expected_random_loci_hits[i])
    if (any(locus >= exon_coordinates$Start & locus <= exon_coordinates$End)) {
      random_exon_count <- random_exon_count + 1
    } else {
      random_intron_count <- random_intron_count + 1 
    }   
    
  }
  
  random_exon_count_summary[y] <- random_exon_count
  random_intron_count_summary[y] <- random_intron_count
  
}

random_exon_count_mean <- mean(unlist(random_exon_count_summary))
random_intron_count_mean <- mean(unlist(random_intron_count_summary))


exon_count_stat <- t.test(unlist(random_exon_count_summary))
intron_count_stat <- t.test(unlist(random_intron_count_summary))


hist(unlist(random_exon_count_summary))
hist(unlist(random_intron_count_summary))


#If SNVs/Indels split as per exon:intron bp ratio 

RAD51B_exon_intron_ratio <- 0.002075005

1570/(1+RAD51B_exon_intron_ratio)
#Expected 1567 to be intronic SNVs/Indels and 3 to be exonic 

RAD51B_exon
RAD51B_intron





