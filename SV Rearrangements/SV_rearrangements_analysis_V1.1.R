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

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/SV Rearrangements/')

WG_SV_RAD51B_Exon_Mut <- as.data.frame(fread('PCAWG_SV_Genome_Mutations_RAD51B_Exon_Donors.tsv', header = TRUE))


WG_SV_RAD51B_Exon_Mut$SV_loci_summary <- paste(WG_SV_RAD51B_Exon_Mut$chr_from, ':', WG_SV_RAD51B_Exon_Mut$chr_from_bkpt, '-', WG_SV_RAD51B_Exon_Mut$chr_to, ':', WG_SV_RAD51B_Exon_Mut$chr_to_bkpt, sep='')
View(WG_SV_RAD51B_Exon_Mut)
nrow(WG_SV_RAD51B_Exon_Mut)
View(table(WG_SV_RAD51B_Exon_Mut$SV_loci_summary))

freq_SV_rearrangements <- WG_SV_RAD51B_Exon_Mut[ WG_SV_RAD51B_Exon_Mut$SV_loci_summary %in% names(table(WG_SV_RAD51B_Exon_Mut$SV_loci_summary))[table(WG_SV_RAD51B_Exon_Mut$SV_loci_summary) >1] , ]


View(freq_SV_rearrangements) #1734
View(WG_SV_RAD51B_Exon_Mut) #4325

length(unique(WG_SV_RAD51B_Exon_Mut$sv_id))

colnames(freq_SV_rearrangements)

View(table(freq_SV_rearrangements$SV_loci_summary))
length(unique(freq_SV_rearrangements$sv_id))
View(table(freq_SV_rearrangements$variant_type))

ofer_SV_calls <- read.csv('merged_1.6.1.csv')


View(ofer_SV_calls)
ofer_SV_calls$SV_loci_summary <- paste(ofer_SV_calls$seqnames, ':', ofer_SV_calls$start, '-', ofer_SV_calls$altchr, ':', ofer_SV_calls$altpos, sep='')
nrow(ofer_SV_calls)
View(table(ofer_SV_calls$SV_loci_summary))

ofer_freq_SV_rearrangements <- ofer_SV_calls[ ofer_SV_calls$SV_loci_summary %in% names(table(ofer_SV_calls$SV_loci_summary))[table(ofer_SV_calls$SV_loci_summary) >1] , ]
