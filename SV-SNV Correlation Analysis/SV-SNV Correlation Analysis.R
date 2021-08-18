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

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/SV-SNV Correlation Analysis/')


SV_calls <- as.data.frame(fread('PCAWG_SV_Genome_Mutations_RAD51B_Exon_Donors.tsv', header = TRUE))

View(SV_calls)

SV_calls$start_locus_summary <-paste(SV_calls$chr_from, SV_calls$chr_from_bkpt,sep = ':')


SV_calls$start_less_50bp <- SV_calls$chr_from_bkpt-50
SV_calls$start_less_50bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_less_50bp,sep = ':')

SV_calls$start_less_100bp <- SV_calls$chr_from_bkpt-100
SV_calls$start_less_100bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_less_100bp,sep = ':')

SV_calls$start_less_250bp <- SV_calls$chr_from_bkpt-250
SV_calls$start_less_250bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_less_250bp,sep = ':')
  
SV_calls$start_less_500bp <- SV_calls$chr_from_bkpt-500
SV_calls$start_less_500bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_less_500bp,sep = ':')
  
SV_calls$start_less_1000bp <- SV_calls$chr_from_bkpt-1000
SV_calls$start_less_1000bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_less_1000bp,sep = ':')
  

SV_calls$start_more_50bp <- SV_calls$chr_from_bkpt+50
SV_calls$start_more_50bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_more_50bp,sep = ':')

SV_calls$start_more_100bp <- SV_calls$chr_from_bkpt+100
SV_calls$start_more_100bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_more_100bp,sep = ':')

SV_calls$start_more_250bp <- SV_calls$chr_from_bkpt+250
SV_calls$start_more_250bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_more_250bp,sep = ':')

SV_calls$start_more_500bp <- SV_calls$chr_from_bkpt+500
SV_calls$start_more_500bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_more_500bp,sep = ':')

SV_calls$start_more_1000bp <- SV_calls$chr_from_bkpt+1000
SV_calls$start_more_1000bp_locus_summary <-paste(SV_calls$chr_from,SV_calls$start_more_1000bp,sep = ':')

  

  

SV_calls$end_locus_summary <- paste(SV_calls$chr_to,SV_calls$chr_to_bkpt, sep=':')


SV_calls$altpos_less_50bp <- SV_calls$chr_to_bkpt-50
SV_calls$altpos_less_50bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_less_50bp,sep = ':')

SV_calls$altpos_less_100bp <- SV_calls$chr_to_bkpt-100
SV_calls$altpos_less_100bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_less_100bp,sep = ':')

SV_calls$altpos_less_250bp <- SV_calls$chr_to_bkpt-250
SV_calls$altpos_less_250bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_less_250bp,sep = ':')

SV_calls$altpos_less_500bp <- SV_calls$chr_to_bkpt-500
SV_calls$altpos_less_500bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_less_500bp,sep = ':')

SV_calls$altpos_less_1000bp <- SV_calls$chr_to_bkpt-1000
SV_calls$altpos_less_1000bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_less_1000bp,sep = ':')



SV_calls$altpos_more_50bp <- SV_calls$chr_to_bkpt+50
SV_calls$altpos_more_50bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_more_50bp,sep = ':')

SV_calls$altpos_more_100bp <- SV_calls$chr_to_bkpt+100
SV_calls$altpos_more_100bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_more_100bp,sep = ':')

SV_calls$altpos_more_250bp <- SV_calls$chr_to_bkpt+250
SV_calls$altpos_more_250bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_more_250bp,sep = ':')

SV_calls$altpos_more_500bp <- SV_calls$chr_to_bkpt+500
SV_calls$altpos_more_500bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_more_500bp,sep = ':')

SV_calls$altpos_more_1000bp <- SV_calls$chr_to_bkpt+1000
SV_calls$altpos_more_1000bp_locus_summary <-paste(SV_calls$chr_to,SV_calls$altpos_more_1000bp,sep = ':')






memory.limit(1e9)



RAD51B_exon_mutation_SNV_Indels <- as.data.frame(fread('PCAWG_SNV_Genome_Mutations_RAD51B_Exon_Donors.tsv', header = TRUE))


RAD51B_exon_mutation_SNV_Indels$locus_summary <- paste(RAD51B_exon_mutation_SNV_Indels$chromosome, RAD51B_exon_mutation_SNV_Indels$chromosome_start, sep=':')
RAD51B_exon_mutation_SNV_Indels$mutation_summary <- paste(RAD51B_exon_mutation_SNV_Indels$mutated_from_allele, RAD51B_exon_mutation_SNV_Indels$mutated_to_allele, sep = '>')




SV_bkpt_start_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_locus_summary)

SV_bkpt_start_less_50bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_less_50bp_locus_summary)
SV_bkpt_start_less_100bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_less_100bp_locus_summary)
SV_bkpt_start_less_250bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_less_250bp_locus_summary)
SV_bkpt_start_less_500bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_less_500bp_locus_summary)
SV_bkpt_start_less_1000bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_less_1000bp_locus_summary)


SV_bkpt_start_more_50bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_more_50bp_locus_summary)
SV_bkpt_start_more_100bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_more_100bp_locus_summary)
SV_bkpt_start_more_250bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_more_250bp_locus_summary)
SV_bkpt_start_more_500bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_more_500bp_locus_summary)
SV_bkpt_start_more_1000bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$start_more_1000bp_locus_summary)






SV_bkpt_end_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$end_locus_summary)


SV_bkpt_end_less_50bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_less_50bp_locus_summary)
SV_bkpt_end_less_100bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_less_100bp_locus_summary)
SV_bkpt_end_less_250bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_less_250bp_locus_summary)
SV_bkpt_end_less_500bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_less_500bp_locus_summary)
SV_bkpt_end_less_1000bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_less_1000bp_locus_summary)


SV_bkpt_end_more_50bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_more_50bp_locus_summary)
SV_bkpt_end_more_100bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_more_100bp_locus_summary)
SV_bkpt_end_more_250bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_more_250bp_locus_summary)
SV_bkpt_end_more_500bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_more_500bp_locus_summary)
SV_bkpt_end_more_1000bp_SNV_Indel <- intersect(RAD51B_exon_mutation_SNV_Indels$locus_summary, SV_calls$altpos_more_1000bp_locus_summary)


SV_bkpt_end_SNV_Indel_summary <- subset(RAD51B_exon_mutation_SNV_Indels, locus_summary %in% SV_bkpt_end_SNV_Indel)
View(SV_bkpt_end_SNV_Indel_summary)

write.csv( c( length(unique(SV_bkpt_end_SNV_Indel_summary$icgc_mutation_id)) , length(unique(SV_bkpt_end_SNV_Indel_summary$icgc_donor_id)) , length(unique(SV_bkpt_end_SNV_Indel_summary$icgc_specimen_id)) , unique(SV_bkpt_end_SNV_Indel_summary$chromosome) , length(unique(SV_bkpt_end_SNV_Indel_summary$icgc_sample_id)) ) , file = paste(deparse(substitute(x)),'_summary_stats','.csv',sep=''))


loci_intersect_summary <- function(x) {
  summary_table <- subset(RAD51B_exon_mutation_SNV_Indels, locus_summary %in% x)
  write.csv(as.data.frame(table(summary_table$mutation_summary)), file = paste(deparse(substitute(x)),'_mutation_from_to','.csv',sep=''))
  write.csv( c( length(unique(summary_table$icgc_mutation_id)) , length(unique(summary_table$icgc_donor_id)) , length(unique(summary_table$icgc_specimen_id)) , unique(summary_table$chromosome) , length(unique(summary_table$icgc_sample_id)) ) , file = paste(deparse(substitute(x)),'_summary_stats','.csv',sep=''))
  write.csv(as.data.frame(table(summary_table$mutation_type)), file = paste(deparse(substitute(x)),'_mutation_type','.csv',sep=''))
  write.csv(as.data.frame(table(summary_table$consequence_type)), file = paste(deparse(substitute(x)),'_consequence_type','.csv',sep=''))
  write.csv(as.data.frame(table(summary_table$gene_affected)), file = paste(deparse(substitute(x)),'_gene_affected','.csv',sep=''))
  write.csv(as.data.frame(table(summary_table$transcript_affected)), file = paste(deparse(substitute(x)),'_transcript_affected','.csv',sep=''))
}


loci_intersect_summary(SV_bkpt_start_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_less_50bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_less_100bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_less_250bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_less_500bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_less_1000bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_more_50bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_more_100bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_more_250bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_more_500bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_start_more_1000bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_less_50bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_less_100bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_less_250bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_less_500bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_less_1000bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_more_50bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_more_100bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_more_250bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_more_500bp_SNV_Indel)
loci_intersect_summary(SV_bkpt_end_more_1000bp_SNV_Indel)

