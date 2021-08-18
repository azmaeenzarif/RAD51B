library(Biostrings)
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

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Intron Comparison/HR Genes')


#HR genes 

RAD51_exon <- readDNAStringSet('exon_Homo_sapiens_RAD51_sequence.fa')
RAD51B_exon <- readDNAStringSet('exon_Homo_sapiens_RAD51B_sequence.fa')
RAD51D_exon <- readDNAStringSet('exon_Homo_sapiens_RAD51D_sequence.fa')
DMC1_exon <- readDNAStringSet('exon_Homo_sapiens_DMC1_sequence.fa')
XRCC2_exon <- readDNAStringSet('exon_Homo_sapiens_XRCC2_sequence.fa')
XRCC3_exon <- readDNAStringSet('exon_Homo_sapiens_XRCC3_sequence.fa')
RAD52_exon <- readDNAStringSet('exon_Homo_sapiens_RAD52_sequence.fa')
RAD54L_exon <- readDNAStringSet('exon_Homo_sapiens_RAD54L_sequence.fa')
RAD54B_exon <- readDNAStringSet('exon_Homo_sapiens_RAD54B_sequence.fa')
BRCA1_exon <- readDNAStringSet('exon_Homo_sapiens_BRCA1_sequence.fa')
RAD50_exon <- readDNAStringSet('exon_Homo_sapiens_RAD50_sequence.fa')
MRE11_exon <- readDNAStringSet('exon_Homo_sapiens_MRE11_sequence.fa')
NBN_exon <- readDNAStringSet('exon_Homo_sapiens_NBN_sequence.fa')
RBBP8_exon <- readDNAStringSet('exon_Homo_sapiens_RBBP8_sequence.fa')
MUS81_exon <- readDNAStringSet('exon_Homo_sapiens_MUS81_sequence.fa')
EME1_exon <- readDNAStringSet('exon_Homo_sapiens_EME1_sequence.fa')
EME2_exon <- readDNAStringSet('exon_Homo_sapiens_EME2_sequence.fa')
GEN1_exon <- readDNAStringSet('exon_Homo_sapiens_GEN1_sequence.fa')


RAD51_intron <- readDNAStringSet('intron_Homo_sapiens_RAD51_sequence.fa')
RAD51B_intron <- readDNAStringSet('intron_Homo_sapiens_RAD51B_sequence.fa')
RAD51D_intron <- readDNAStringSet('intron_Homo_sapiens_RAD51D_sequence.fa')
DMC1_intron <- readDNAStringSet('intron_Homo_sapiens_DMC1_sequence.fa')
XRCC2_intron <- readDNAStringSet('intron_Homo_sapiens_XRCC2_sequence.fa')
XRCC3_intron <- readDNAStringSet('intron_Homo_sapiens_XRCC3_sequence.fa')
RAD52_intron <- readDNAStringSet('intron_Homo_sapiens_RAD52_sequence.fa')
RAD54L_intron <- readDNAStringSet('intron_Homo_sapiens_RAD54L_sequence.fa')
RAD54B_intron <- readDNAStringSet('intron_Homo_sapiens_RAD54B_sequence.fa')
BRCA1_intron <- readDNAStringSet('intron_Homo_sapiens_BRCA1_sequence.fa')
RAD50_intron <- readDNAStringSet('intron_Homo_sapiens_RAD50_sequence.fa')
MRE11_intron <- readDNAStringSet('intron_Homo_sapiens_MRE11_sequence.fa')
NBN_intron <- readDNAStringSet('intron_Homo_sapiens_NBN_sequence.fa')
RBBP8_intron <- readDNAStringSet('intron_Homo_sapiens_RBBP8_sequence.fa')
MUS81_intron <- readDNAStringSet('intron_Homo_sapiens_MUS81_sequence.fa')
EME1_intron <- readDNAStringSet('intron_Homo_sapiens_EME1_sequence.fa')
EME2_intron <- readDNAStringSet('intron_Homo_sapiens_EME2_sequence.fa')
GEN1_intron <- readDNAStringSet('intron_Homo_sapiens_GEN1_sequence.fa')


HR_genes <- list("RAD51","RAD51B","RAD51D","DMC1","XRCC2","XRCC3","RAD52","RAD54L","RAD54B","BRCA1","RAD50","MRE11","NBN","RBBP8","MUS81","EME1","EME2","GEN1")

intron_exon_comparison_table <- data.frame()




for (i in 1:length(HR_genes)) {
  gene_exon <- as.name(gsub(' ', '', paste(HR_genes[i],'_exon')))
  gene_intron <- as.name(gsub(' ', '', paste(HR_genes[i],'_intron')))
  intron_exon_comparison_table[i,1] <- toString(HR_genes[i])
  intron_exon_comparison_table[i,2] <- ((nrow(as.data.frame(unlist(eval(gene_exon)))))/(nrow(as.data.frame(unlist(eval(gene_intron))))))
}


colnames(intron_exon_comparison_table) <- c('HR Gene', 'Exon:Intron Ratio per intronic bp' )



View(intron_exon_comparison_table)


intron_exon_comparison_table <- intron_exon_comparison_table[order(intron_exon_comparison_table$`Exon:Intron Ratio per intronic bp`, decreasing= TRUE),]


barplot(intron_exon_comparison_table$`Exon:Intron Ratio per intronic bp`, names.arg = intron_exon_comparison_table$`HR Gene`, las=2, cex.names = 0.6, ylab = 'Exon:Intron Ratio (per bp)', xlab = 'HR Gene')



#Common Cancer Genes

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Intron Comparison/Cancer Genes')


ABL1_exon <- readDNAStringSet('exon_Homo_sapiens_ABL1_sequence.fa')
AFF4_exon <- readDNAStringSet('exon_Homo_sapiens_AFF4_sequence.fa')
AKAP13_exon <- readDNAStringSet('exon_Homo_sapiens_AKAP13_sequence.fa')
AKT2_exon <- readDNAStringSet('exon_Homo_sapiens_AKT2_sequence.fa')
ALK_exon <- readDNAStringSet('exon_Homo_sapiens_ALK_sequence.fa')
APC_exon <- readDNAStringSet('exon_Homo_sapiens_APC_sequence.fa')
AXL_exon <- readDNAStringSet('exon_Homo_sapiens_AXL_sequence.fa')
BCL2_exon <- readDNAStringSet('exon_Homo_sapiens_BCL2_sequence.fa')
#BRCA1_exon <- readDNAStringSet('exon_Homo_sapiens_BRCA1_sequence.fa')
BRCA2_exon <- readDNAStringSet('exon_Homo_sapiens_BRCA2_sequence.fa')
CCND1_exon <- readDNAStringSet('exon_Homo_sapiens_CCND1_sequence.fa')
CDKN2A_exon <- readDNAStringSet('exon_Homo_sapiens_CDKN2A_sequence.fa')
CSF1R_exon <- readDNAStringSet('exon_Homo_sapiens_CSF1R_sequence.fa')
DCC_exon <- readDNAStringSet('exon_Homo_sapiens_DCC_sequence.fa')
EGFR_exon <- readDNAStringSet('exon_Homo_sapiens_EGFR_sequence.fa')
ERBB2_exon <- readDNAStringSet('exon_Homo_sapiens_ERBB2_sequence.fa')
ETS1_exon <- readDNAStringSet('exon_Homo_sapiens_ETS1_sequence.fa')
FES_exon <- readDNAStringSet('exon_Homo_sapiens_FES_sequence.fa')
FGF3_exon <- readDNAStringSet('exon_Homo_sapiens_FGF3_sequence.fa')
FGF4_exon <- readDNAStringSet('exon_Homo_sapiens_FGF4_sequence.fa')
FOS_exon <- readDNAStringSet('exon_Homo_sapiens_FOS_sequence.fa')
GLI1_exon <- readDNAStringSet('exon_Homo_sapiens_GLI1_sequence.fa')
GNAS_exon <- readDNAStringSet('exon_Homo_sapiens_GNAS_sequence.fa')
HRAS_exon <- readDNAStringSet('exon_Homo_sapiens_HRAS_sequence.fa')
IL3_exon <- readDNAStringSet('exon_Homo_sapiens_IL3_sequence.fa')
JUN_exon <- readDNAStringSet('exon_Homo_sapiens_JUN_sequence.fa')
KIT_exon <- readDNAStringSet('exon_Homo_sapiens_KIT_sequence.fa')
KRAS_exon <- readDNAStringSet('exon_Homo_sapiens_KRAS_sequence.fa')
LCK_exon <- readDNAStringSet('exon_Homo_sapiens_LCK_sequence.fa')
LMO1_exon <- readDNAStringSet('exon_Homo_sapiens_LMO1_sequence.fa')
LMO2_exon <- readDNAStringSet('exon_Homo_sapiens_LMO2_sequence.fa')
LYL1_exon <- readDNAStringSet('exon_Homo_sapiens_LYL1_sequence.fa')
MAS1_exon <- readDNAStringSet('exon_Homo_sapiens_MAS1_sequence.fa')
MCF2L_exon <- readDNAStringSet('exon_Homo_sapiens_MCF2L_sequence.fa')
MDM2_exon <- readDNAStringSet('exon_Homo_sapiens_MDM2_sequence.fa')
MEN1_exon <- readDNAStringSet('exon_Homo_sapiens_MEN1_sequence.fa')
#MOS1_exon <- readDNAStringSet('exon_Homo_sapiens_MOS1_sequence.fa')
MYB_exon <- readDNAStringSet('exon_Homo_sapiens_MYB_sequence.fa')
MYC_exon <- readDNAStringSet('exon_Homo_sapiens_MYC_sequence.fa')
MYCL_exon <- readDNAStringSet('exon_Homo_sapiens_MYCL_sequence.fa')
MYCN_exon <- readDNAStringSet('exon_Homo_sapiens_MYCN_sequence.fa')
NF1_exon <- readDNAStringSet('exon_Homo_sapiens_NF1_sequence.fa')
NF2_exon <- readDNAStringSet('exon_Homo_sapiens_NF2_sequence.fa')
NKFB2_exon <- readDNAStringSet('exon_Homo_sapiens_NFKB2_sequence.fa')
NOTCH1_exon <- readDNAStringSet('exon_Homo_sapiens_NOTCH1_sequence.fa')
NRAS_exon <- readDNAStringSet('exon_Homo_sapiens_NRAS_sequence.fa')
NTRK1_exon <- readDNAStringSet('exon_Homo_sapiens_NTRK1_sequence.fa')
PAX5_exon <- readDNAStringSet('exon_Homo_sapiens_PAX5_sequence.fa')
PDGFB_exon <- readDNAStringSet('exon_Homo_sapiens_PDGFB_sequence.fa')
PIM1_exon <- readDNAStringSet('exon_Homo_sapiens_PIM1_sequence.fa')
PTEN_exon <- readDNAStringSet('exon_Homo_sapiens_PTEN_sequence.fa')
RAF1_exon <- readDNAStringSet('exon_Homo_sapiens_RAF1_sequence.fa')
RB1_exon <- readDNAStringSet('exon_Homo_sapiens_RB1_sequence.fa')
#RET1_exon <- readDNAStringSet('exon_Homo_sapiens_RET1_sequence.fa')
ROS1_exon <- readDNAStringSet('exon_Homo_sapiens_ROS1_sequence.fa')
RUNX1_exon <- readDNAStringSet('exon_Homo_sapiens_RUNX1_sequence.fa')
SKI_exon <- readDNAStringSet('exon_Homo_sapiens_SKI_sequence.fa')
SMAD2_exon <- readDNAStringSet('exon_Homo_sapiens_SMAD2_sequence.fa')
SMAD4_exon <- readDNAStringSet('exon_Homo_sapiens_SMAD4_sequence.fa')
SRC_exon <- readDNAStringSet('exon_Homo_sapiens_SRC_sequence.fa')
TAL1_exon <- readDNAStringSet('exon_Homo_sapiens_TAL1_sequence.fa')
TAL2_exon <- readDNAStringSet('exon_Homo_sapiens_TAL2_sequence.fa')
TIAM1_exon <- readDNAStringSet('exon_Homo_sapiens_TIAM1_sequence.fa')
TLX1_exon <- readDNAStringSet('exon_Homo_sapiens_TLX1_sequence.fa')
TP53_exon <- readDNAStringSet('exon_Homo_sapiens_TP53_sequence.fa')
TSC2_exon <- readDNAStringSet('exon_Homo_sapiens_TSC2_sequence.fa')
VHL_exon <- readDNAStringSet('exon_Homo_sapiens_VHL_sequence.fa')
WRN_exon <- readDNAStringSet('exon_Homo_sapiens_WRN_sequence.fa')
WT1_exon <- readDNAStringSet('exon_Homo_sapiens_WT1_sequence.fa')

RAD51B_exon <- readDNAStringSet('exon_Homo_sapiens_RAD51B_sequence.fa')



ABL1_intron <- readDNAStringSet('intron_Homo_sapiens_ABL1_sequence.fa')
AFF4_intron <- readDNAStringSet('intron_Homo_sapiens_AFF4_sequence.fa')
AKAP13_intron <- readDNAStringSet('intron_Homo_sapiens_AKAP13_sequence.fa')
AKT2_intron <- readDNAStringSet('intron_Homo_sapiens_AKT2_sequence.fa')
ALK_intron <- readDNAStringSet('intron_Homo_sapiens_ALK_sequence.fa')
APC_intron <- readDNAStringSet('intron_Homo_sapiens_APC_sequence.fa')
AXL_intron <- readDNAStringSet('intron_Homo_sapiens_AXL_sequence.fa')
BCL2_intron <- readDNAStringSet('intron_Homo_sapiens_BCL2_sequence.fa')
#BRCA1_intron <- readDNAStringSet('intron_Homo_sapiens_BRCA1_sequence.fa')
BRCA2_intron <- readDNAStringSet('intron_Homo_sapiens_BRCA2_sequence.fa')
CCND1_intron <- readDNAStringSet('intron_Homo_sapiens_CCND1_sequence.fa')
CDKN2A_intron <- readDNAStringSet('intron_Homo_sapiens_CDKN2A_sequence.fa')
CSF1R_intron <- readDNAStringSet('intron_Homo_sapiens_CSF1R_sequence.fa')
DCC_intron <- readDNAStringSet('intron_Homo_sapiens_DCC_sequence.fa')
EGFR_intron <- readDNAStringSet('intron_Homo_sapiens_EGFR_sequence.fa')
ERBB2_intron <- readDNAStringSet('intron_Homo_sapiens_ERBB2_sequence.fa')
ETS1_intron <- readDNAStringSet('intron_Homo_sapiens_ETS1_sequence.fa')
FES_intron <- readDNAStringSet('intron_Homo_sapiens_FES_sequence.fa')
FGF3_intron <- readDNAStringSet('intron_Homo_sapiens_FGF3_sequence.fa')
FGF4_intron <- readDNAStringSet('intron_Homo_sapiens_FGF4_sequence.fa')
FOS_intron <- readDNAStringSet('intron_Homo_sapiens_FOS_sequence.fa')
GLI1_intron <- readDNAStringSet('intron_Homo_sapiens_GLI1_sequence.fa')
GNAS_intron <- readDNAStringSet('intron_Homo_sapiens_GNAS_sequence.fa')
HRAS_intron <- readDNAStringSet('intron_Homo_sapiens_HRAS_sequence.fa')
IL3_intron <- readDNAStringSet('intron_Homo_sapiens_IL3_sequence.fa')
JUN_intron <- readDNAStringSet('intron_Homo_sapiens_JUN_sequence.fa')
KIT_intron <- readDNAStringSet('intron_Homo_sapiens_KIT_sequence.fa')
KRAS_intron <- readDNAStringSet('intron_Homo_sapiens_KRAS_sequence.fa')
LCK_intron <- readDNAStringSet('intron_Homo_sapiens_LCK_sequence.fa')
LMO1_intron <- readDNAStringSet('intron_Homo_sapiens_LMO1_sequence.fa')
LMO2_intron <- readDNAStringSet('intron_Homo_sapiens_LMO2_sequence.fa')
LYL1_intron <- readDNAStringSet('intron_Homo_sapiens_LYL1_sequence.fa')
MAS1_intron <- readDNAStringSet('intron_Homo_sapiens_MAS1_sequence.fa')
MCF2L_intron <- readDNAStringSet('intron_Homo_sapiens_MCF2L_sequence.fa')
MDM2_intron <- readDNAStringSet('intron_Homo_sapiens_MDM2_sequence.fa')
MEN1_intron <- readDNAStringSet('intron_Homo_sapiens_MEN1_sequence.fa')
#MOS1_intron <- readDNAStringSet('intron_Homo_sapiens_MOS1_sequence.fa')
MYB_intron <- readDNAStringSet('intron_Homo_sapiens_MYB_sequence.fa')
MYC_intron <- readDNAStringSet('intron_Homo_sapiens_MYC_sequence.fa')
MYCL_intron <- readDNAStringSet('intron_Homo_sapiens_MYCL_sequence.fa')
MYCN_intron <- readDNAStringSet('intron_Homo_sapiens_MYCN_sequence.fa')
NF1_intron <- readDNAStringSet('intron_Homo_sapiens_NF1_sequence.fa')
NF2_intron <- readDNAStringSet('intron_Homo_sapiens_NF2_sequence.fa')
NKFB2_intron <- readDNAStringSet('intron_Homo_sapiens_NFKB2_sequence.fa')
NOTCH1_intron <- readDNAStringSet('intron_Homo_sapiens_NOTCH1_sequence.fa')
NRAS_intron <- readDNAStringSet('intron_Homo_sapiens_NRAS_sequence.fa')
NTRK1_intron <- readDNAStringSet('intron_Homo_sapiens_NTRK1_sequence.fa')
PAX5_intron <- readDNAStringSet('intron_Homo_sapiens_PAX5_sequence.fa')
PDGFB_intron <- readDNAStringSet('intron_Homo_sapiens_PDGFB_sequence.fa')
PIM1_intron <- readDNAStringSet('intron_Homo_sapiens_PIM1_sequence.fa')
PTEN_intron <- readDNAStringSet('intron_Homo_sapiens_PTEN_sequence.fa')
RAF1_intron <- readDNAStringSet('intron_Homo_sapiens_RAF1_sequence.fa')
RB1_intron <- readDNAStringSet('intron_Homo_sapiens_RB1_sequence.fa')
#RET1_intron <- readDNAStringSet('intron_Homo_sapiens_RET1_sequence.fa')
ROS1_intron <- readDNAStringSet('intron_Homo_sapiens_ROS1_sequence.fa')
RUNX1_intron <- readDNAStringSet('intron_Homo_sapiens_RUNX1_sequence.fa')
SKI_intron <- readDNAStringSet('intron_Homo_sapiens_SKI_sequence.fa')
SMAD2_intron <- readDNAStringSet('intron_Homo_sapiens_SMAD2_sequence.fa')
SMAD4_intron <- readDNAStringSet('intron_Homo_sapiens_SMAD4_sequence.fa')
SRC_intron <- readDNAStringSet('intron_Homo_sapiens_SRC_sequence.fa')
TAL1_intron <- readDNAStringSet('intron_Homo_sapiens_TAL1_sequence.fa')
TAL2_intron <- readDNAStringSet('intron_Homo_sapiens_TAL2_sequence.fa')
TIAM1_intron <- readDNAStringSet('intron_Homo_sapiens_TIAM1_sequence.fa')
TLX1_intron <- readDNAStringSet('intron_Homo_sapiens_TLX1_sequence.fa')
TP53_intron <- readDNAStringSet('intron_Homo_sapiens_TP53_sequence.fa')
TSC2_intron <- readDNAStringSet('intron_Homo_sapiens_TSC2_sequence.fa')
VHL_intron <- readDNAStringSet('intron_Homo_sapiens_VHL_sequence.fa')
WRN_intron <- readDNAStringSet('intron_Homo_sapiens_WRN_sequence.fa')
WT1_intron <- readDNAStringSet('intron_Homo_sapiens_WT1_sequence.fa')

RAD51B_intron <- readDNAStringSet('intron_Homo_sapiens_RAD51B_sequence.fa')



Cancer_genes <- list('ABL1','AFF4','AKAP13','AKT2','ALK','APC','AXL','BCL2','BRCA2','CCND1','CSF1R','DCC','EGFR','ERBB2','ETS1','FES','FGF3','FGF4','FOS','GLI1','GNAS','HRAS','IL3','JUN','KIT','KRAS','LCK','LMO1','LMO2','LYL1','MAS1','MCF2L','MDM2','MEN1','MYB','MYC','MYCL','MYCN','NF1','NF2','NKFB2','NOTCH1','NRAS','NTRK1','PAX5','PDGFB','PIM1','PTEN','RAF1','RB1','ROS1','RUNX1','SKI','SMAD2','SMAD4','SRC','TAL1','TIAM1','TLX1','TP53','TSC2','VHL','WRN','WT1', 'RAD51B')
  
  
cancer_genes_intron_exon_comparison_table <- data.frame()




for (i in 1:length(Cancer_genes)) {
  gene_exon <- as.name(gsub(' ', '', paste(Cancer_genes[i],'_exon')))
  gene_intron <- as.name(gsub(' ', '', paste(Cancer_genes[i],'_intron')))
  cancer_genes_intron_exon_comparison_table[i,1] <- toString(Cancer_genes[i])
  cancer_genes_intron_exon_comparison_table[i,2] <- ((nrow(as.data.frame(unlist(eval(gene_exon)))))/(nrow(as.data.frame(unlist(eval(gene_intron))))))
}


colnames(cancer_genes_intron_exon_comparison_table) <- c('Cancer Gene', 'Exon:Intron Ratio per intronic bp' )

cancer_genes_intron_exon_comparison_table <- cancer_genes_intron_exon_comparison_table[order(cancer_genes_intron_exon_comparison_table$`Exon:Intron Ratio per intronic bp`, decreasing= TRUE),]


barplot(cancer_genes_intron_exon_comparison_table$`Exon:Intron Ratio per intronic bp`, names.arg = cancer_genes_intron_exon_comparison_table$`Cancer Gene`, las=2, cex.names = 0.5, ylab = 'Exon:Intron Ratio (per bp)', xlab = 'Known Cancer Gene')

