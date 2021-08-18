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
library(readr)

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Question 2/')


exon_mutations <- read.csv('RAD51B_exon_mutations_final.csv')
exon_mutations$no_of_donors <- substring(exon_mutations$Donors.affected, 3)
View(exon_mutations)


mut_ids <- data.frame(exon_mutations$Mutation.ID)
#View(mut_ids)

exon_mutations_donor_ids <- list()


SGV_count <- 0
expArray_count <- 0
expSeq_count <- 0



get_donors <- GET(gsub(' ', '', paste('https://dcc.icgc.org/api/v1/mutations/MU2395912/donors')))
get_donors_text <- content(get_donors, 'text', encoding = 'UTF-8')
get_donors_json <- fromJSON(get_donors_text, flatten=TRUE)
get_donors_json_df <- data.frame(get_donors_json = unlist(get_donors_json))
View(get_donors_json_df[1])




for (i in 1:nrow(mut_ids)) {
  mutation <- mut_ids[i,]
  get_donors <- GET(gsub(' ', '', paste('https://dcc.icgc.org/api/v1/mutations/',mutation, '/donors')))
  get_donors_text <- content(get_donors, 'text', encoding = 'UTF-8')
  get_donors_json <- fromJSON(get_donors_text, flatten=TRUE)
  get_donors_json_df <- data.frame(get_donors_json = unlist(get_donors_json))
  
  
  if (exon_mutations$no_of_donors[i] == 1) {
    exon_mutations_donor_ids <- c(exon_mutations_donor_ids, get_donors_json_df[1,])
    if (get_donors_json_df['hits.sgvExists',] == TRUE) {
      SGV_count <- SGV_count+1
    }
    
    if (get_donors_json_df['hits.expSeqExists',] == TRUE) {
      expSeq_count <- expSeq_count+1
    }
    
    if (get_donors_json_df['hits.expArrayExists',] == TRUE) {
      expArray_count <- expArray_count+1
    }
  
    
    
  }
  
  else {
    if (exon_mutations$no_of_donors[i] == 2) {
      exon_mutations_donor_ids <- c(exon_mutations_donor_ids, get_donors_json_df[1:2,])
      
      
      if (get_donors_json_df['hits.sgvExists1',] == TRUE) {
        SGV_count <- SGV_count+1
      }
      
      if (get_donors_json_df['hits.expSeqExists1',] == TRUE) {
        expSeq_count <- expSeq_count+1
      }
      
      if (get_donors_json_df['hits.expArrayExists1',] == TRUE) {
        expArray_count <- expArray_count+1
      }
      
      
      if (get_donors_json_df['hits.sgvExists2',] == TRUE) {
        SGV_count <- SGV_count+1
      }
      
      if (get_donors_json_df['hits.expSeqExists2',] == TRUE) {
        expSeq_count <- expSeq_count+1
      }
      
      if (get_donors_json_df['hits.expArrayExists2',] == TRUE) {
        expArray_count <- expArray_count+1
      }
      
      
      
      
    }
    
    else {
  
      if (exon_mutations$no_of_donors[i] == 3) {
        exon_mutations_donor_ids <- c(exon_mutations_donor_ids, get_donors_json_df[1:3,])
        if (get_donors_json_df['hits.sgvExists1',] == TRUE) {
          SGV_count <- SGV_count+1
        }
        
        if (get_donors_json_df['hits.expSeqExists1',] == TRUE) {
          expSeq_count <- expSeq_count+1
        }
        
        if (get_donors_json_df['hits.expArrayExists1',] == TRUE) {
          expArray_count <- expArray_count+1
        }
        
        if (get_donors_json_df['hits.sgvExists2',] == TRUE) {
          SGV_count <- SGV_count+1
        }
        
        if (get_donors_json_df['hits.expSeqExists2',] == TRUE) {
          expSeq_count <- expSeq_count+1
        }
        
        if (get_donors_json_df['hits.expArrayExists2',] == TRUE) {
          expArray_count <- expArray_count+1
        }
        
        
        if (get_donors_json_df['hits.sgvExists3',] == TRUE) {
          SGV_count <- SGV_count+1
        }
        
        if (get_donors_json_df['hits.expSeqExists3',] == TRUE) {
          expSeq_count <- expSeq_count+1
        }
        
        if (get_donors_json_df['hits.expArrayExists3',] == TRUE) {
          expArray_count <- expArray_count+1
        }
        
      }
    }
  }
}



SGV_count
expArray_count
expSeq_count



exon_mutations$no_of_donors <- as.numeric(exon_mutations$no_of_donors)

exon_mutations_donor_ids <- as.data.frame(exon_mutations_donor_ids)


write.csv(exon_mutations_donor_ids, 'exon_mutations_donor_ids.csv')




memory.limit(size=1000000000)


WG_CNV_RAD51B_Exon_Mut <- as.data.frame(fread('PCAWG_CNV_Genome_Mutations_RAD51B_Exon_Donors.tsv', header = TRUE))
WG_SNV_RAD51B_Exon_Mut <- as.data.frame(fread('PCAWG_SNV_Genome_Mutations_RAD51B_Exon_Donors.tsv', header = TRUE))
WG_SV_RAD51B_Exon_Mut <- as.data.frame(fread('PCAWG_SV_Genome_Mutations_RAD51B_Exon_Donors.tsv', header = TRUE))



View(WG_CNV_RAD51B_Exon_Mut)
View(WG_SNV_RAD51B_Exon_Mut)
View(WG_SV_RAD51B_Exon_Mut)





#SNV Analysis

SNV_mutation_type_freq_table <- table(WG_SNV_RAD51B_Exon_Mut$mutation_type)
View(SNV_mutation_type_freq_table)
SNV_mutations_per_donor_freq_table <- table(WG_SNV_RAD51B_Exon_Mut$icgc_donor_id)

mean_SNV_mutations_per_donor <- mean(SNV_mutations_per_donor_freq_table)
mean_SNV_mutations_per_donor

median_SNV_mutations_per_donor <- median(SNV_mutations_per_donor_freq_table)
median_SNV_mutations_per_donor

SNV_mutations_chrs_affected <- table(WG_SNV_RAD51B_Exon_Mut$chromosome)
View(SNV_mutations_chrs_affected)



SNV_mutation_total_read_count <- table(WG_SNV_RAD51B_Exon_Mut$total_read_count)
mean_SNV_mutation_total_read_count <- mean(SNV_mutation_total_read_count)
median_SNV_mutation_total_read_count <- median(SNV_mutation_total_read_count)


SNV_mutation_mutant_allele_read_count <- table(WG_SNV_RAD51B_Exon_Mut$mutant_allele_read_count)
mean_SNV_mutation_mutant_allele_read_count <- mean(SNV_mutation_mutant_allele_read_count)
median_SNV_mutation_mutant_allele_read_count <- median(SNV_mutation_mutant_allele_read_count)
median_SNV_mutation_mutant_allele_read_count

SNV_mutation_consequency_type_freq_table <- table(WG_SNV_RAD51B_Exon_Mut$consequence_type)
View(SNV_mutation_consequency_type_freq_table)
sum(SNV_mutation_consequency_type_freq_table)

SNV_mutation_gene_affected_freq_table <- as.data.frame(table(WG_SNV_RAD51B_Exon_Mut$gene_affected))
View(SNV_mutation_gene_affected_freq_table)
sum(SNV_mutation_transcript_affected_freq_table)

SNV_mutation_transcript_affected_freq_table <- table(WG_SNV_RAD51B_Exon_Mut$transcript_affected)
View(SNV_mutation_transcript_affected_freq_table)

#SV Analysis

SV_variant_type_freq_table <- table(WG_SV_RAD51B_Exon_Mut$variant_type)
View(SV_variant_type_freq_table)
sum(SV_variant_type_freq_table)

SV_mutations_per_donor_freq_table <- table(WG_SV_RAD51B_Exon_Mut$icgc_donor_id)
mean_SV_mutations_per_donor <- mean(SV_mutations_per_donor_freq_table)
mean_SV_mutations_per_donor
median_SV_mutations_per_donor <- median(SV_mutations_per_donor_freq_table)
median_SV_mutations_per_donor



#CNV Analysis


CNV_mutation_type_freq_table <- table(WG_CNV_RAD51B_Exon_Mut$mutation_type)
View(CNV_mutation_type_freq_table)
sum(CNV_mutation_type_freq_table)

CNV_mutations_per_donor_freq_table <- table(WG_CNV_RAD51B_Exon_Mut$icgc_donor_id)
mean_CNV_mutations_per_donor <- mean(CNV_mutations_per_donor_freq_table)
median_CNV_mutations_per_donor <- median(CNV_mutations_per_donor_freq_table)
median_CNV_mutations_per_donor


CNV_mutations_chrs_affected <- table(WG_CNV_RAD51B_Exon_Mut$chromosome)

CNV_mutations_copy_number_change <- table(WG_CNV_RAD51B_Exon_Mut$copy_number)
View(CNV_mutations_copy_number_change)

CNV_segment_mean_freq_table <- table(sign(WG_CNV_RAD51B_Exon_Mut$segment_mean))
CNV_segment_mean_freq_table




RAD51B_exon_mutation_list <- list.files("RAD51B_Mutations/", pattern = '\\.tsv$', full.names = TRUE)
all_RAD51B_mutations <- lapply(RAD51B_exon_mutation_list, function(x){
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})


all_RAD51B_mutations <- bind_rows(all_RAD51B_mutations)
View(all_RAD51B_mutations)
View(table(all_RAD51B_mutations$Type))


all_RAD51B_mutations$Mutation_Nature <- substring(all_RAD51B_mutations$Genomic.DNA.Change, 17)

View(all_RAD51B_mutations$Mutation_Nature)
View(table(all_RAD51B_mutations$Mutation_Nature))


all_RAD51B_mutations$Mutation_Nature_Frequency <- table(all_RAD51B_mutations$Mutation_Nature)[all_RAD51B_mutations$Mutation_Nature]

all_RAD51B_mutations$Mutation_Nature_Rank <- rank(all_RAD51B_mutations$Mutation_Nature_Frequency, ties.method='min')



###Recurring SNVs/Indels in RAD51B

RAD51B_SNVs <- as.data.frame(subset(WG_SNV_RAD51B_Exon_Mut, chromosome == '14' & chromosome_start >= '67819779' & chromosome_end <= '68730218'))
length(unique(RAD51B_SNVs$icgc_donor_id))
nrow(RAD51B_SNVs)


RAD51B_SNV_insertions <- subset(RAD51B_SNVs, mutation_type == 'insertion of <=200bp')

RAD51B_SNV_deletions <- subset(RAD51B_SNVs, mutation_type == 'deletion of <=200bp')

RAD51B_single_base_substitutions <- subset(RAD51B_SNVs, mutation_type == 'single base substitution')

RAD51B_multiple_base_substitutions <- subset(RAD51B_SNVs, mutation_type == 'multiple base substitution (>=2bp and <=200bp)')


#SNV insertions

View(RAD51B_SNV_insertions)


RAD51B_SNV_insertions %>% summarise_all(n_distinct)

length(unique(RAD51B_SNV_insertions$chromosome_start)) #7 unique start loci for 94 insertions (same gene different transcripts?)
unique(RAD51B_SNV_insertions$chromosome_start)

length(unique(RAD51B_SNV_insertions$chromosome_end)) #7 unique end loci for 94 insertions
unique(RAD51B_SNV_insertions$chromosome_end)

unique(RAD51B_SNV_insertions$mutated_to_allele) #3 - "T"    "TTTT" "GG"  

recurring_insertions <-aggregate(numeric(nrow(RAD51B_SNV_insertions)), RAD51B_SNV_insertions[c('chromosome_start', 'mutated_to_allele')], length)
View(recurring_insertions)

RAD51B_SNV_insertions

#SNV deletions

View(RAD51B_SNV_deletions)

RAD51B_SNV_deletions %>% summarise_all(n_distinct)


length(unique(RAD51B_SNV_deletions$chromosome_start)) #17 unique loci 
unique(RAD51B_SNV_deletions$chromosome_start)

length(unique(RAD51B_SNV_deletions$chromosome_end)) #17 unique loci
unique(RAD51B_SNV_deletions$chromosome_end)

unique(RAD51B_SNV_deletions$mutated_from_allele) #9 - "CT"          "T"           "AGA"         "A"           "ATG"         "AG"          "C"           "GTGTAGAGAGG" "TT"

recurring_deletions <-aggregate(numeric(nrow(RAD51B_SNV_deletions)), RAD51B_SNV_deletions[c('chromosome_start', 'mutated_from_allele')], length)
View(recurring_deletions)


#SNV single base substitutions

View(RAD51B_single_base_substitutions)

RAD51B_single_base_substitutions %>% summarise_all(n_distinct)


length(unique(RAD51B_single_base_substitutions$chromosome_start)) #225 unique loci
unique(RAD51B_single_base_substitutions$chromosome_start)

length(unique(RAD51B_single_base_substitutions$chromosome_end))
unique(RAD51B_single_base_substitutions$chromosome_end)



for (i in 1:nrow(RAD51B_single_base_substitutions)) {
  RAD51B_single_base_substitutions$summary_substitution[i] <- paste(toString(RAD51B_single_base_substitutions$mutated_from_allele[i]), '>', toString(RAD51B_single_base_substitutions$mutated_to_allele[i]))
}


recurring_single_base <- aggregate(numeric(nrow(RAD51B_single_base_substitutions)), RAD51B_single_base_substitutions[c('chromosome_start', 'summary_substitution')], length)
View(recurring_single_base)


#SNV multiple base substitutions

View(RAD51B_multiple_base_substitutions)

RAD51B_multiple_base_substitutions %>% summarise_all(n_distinct)


length(unique(RAD51B_multiple_base_substitutions$chromosome_start)) #4 unique loci
unique(RAD51B_multiple_base_substitutions$chromosome_start)

length(unique(RAD51B_multiple_base_substitutions$chromosome_end)) 
unique(RAD51B_multiple_base_substitutions$chromosome_end)

for (i in 1:nrow(RAD51B_multiple_base_substitutions)) {
  RAD51B_multiple_base_substitutions$summary_substitution[i] <- paste(toString(RAD51B_multiple_base_substitutions$mutated_from_allele[i]), '>', toString(RAD51B_multiple_base_substitutions$mutated_to_allele[i]))
}


recurring_multiple_base <- aggregate(numeric(nrow(RAD51B_multiple_base_substitutions)), RAD51B_multiple_base_substitutions[c('chromosome_start', 'summary_substitution')], length)
View(recurring_multiple_base)






