library(reshape2)
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

setwd('C:/Users/Azmaeen Zarif/Documents/Research/Summer Research Project 2021/Code/Validation Analysis/')

#Loads all DNA Repair [Pathway Ref: R-HSA-73894 Pathway] mutations (in 4 files) from PCAWG dataset
dna_repair_mutations_list <- list.files("DNA_Repair_Genes/", pattern = '\\.tsv$', full.names = TRUE)


all_mutations <- lapply(dna_repair_mutations_list, function(x){
  read.table(file = x,
             sep = '\t',
             fill=TRUE,
             quote ='',
             header = TRUE)
})




#Combines all files into one data frame
all_dna_repair_mutations <- bind_rows(all_mutations)


all_dna_repair_mutations$Affected.Donors <- as.integer(sub('/.*','',all_dna_repair_mutations$Affected.Donors))


for (i in 1:nrow(all_dna_repair_mutations)) {
  all_dna_repair_mutations$percentage_donors_affected[i] <- ((all_dna_repair_mutations$Affected.Donors[i]/2558))
}


###Percentage donor mutations
ggplot(all_dna_repair_mutations, aes(x=Symbol, y=percentage_donors_affected))+
  geom_bar(stat='identity')+
  theme(axis.text.x=element_text(angle =90, vjust = 0.5, hjust = 1))

ggsave('Percentage_dna_repair_mutations_donor.pdf', width = 30, height = 10)

View(all_dna_repair_mutations)


###Absolute number of mutations
ggplot(all_dna_repair_mutations, aes(x=Symbol, y=Mutations))+
  geom_bar(stat='identity')+
  theme(axis.text.x=element_text(angle =90, vjust = 0.5, hjust = 1))
ggsave('Absolute_dna_repair_mutations.pdf', width = 30, height = 10)


###SV Burden Analysis


RAD51B_Not_Mutated_CNV <- data.frame(fread('RAD51B_Not_Mutated_CNV.gz'))
RAD51B_Not_Mutated_Somatic_Mutations <- data.frame(fread('RAD51B_Not_Mutated_Somatic_Mutations.gz'))
RAD51B_Not_Mutated_SV <- data.frame(fread('RAD51B_Not_Mutated_SV.gz'))


RAD51B_Mutated_CNV <- data.frame(fread('RAD51B_Mutated_CNV.gz'))
RAD51B_Mutated_SV <- data.frame(fread('RAD51B_Mutated_SV.gz'))


View(RAD51B_Mutated_SV)

RAD51B_Not_Mutated_SV$icgc_sample_id <- substring(RAD51B_Not_Mutated_SV$icgc_sample_id, 3)
RAD51B_Mutated_SV$icgc_sample_id <- substring(RAD51B_Mutated_SV$icgc_sample_id, 3)

RAD51B_Not_Mutated_SV[order(RAD51B_Not_Mutated_SV$icgc_sample_id),]
RAD51B_Mutated_SV[order(RAD51B_Mutated_SV$icgc_sample_id),]



RAD51B_Not_Mutated_SV_per_Sample <- data.table(RAD51B_Not_Mutated_SV)[,.N, by=RAD51B_Not_Mutated_SV$icgc_sample_id]
setnames(RAD51B_Not_Mutated_SV_per_Sample, 'N', 'RAD51B WT SV Mutation Frequency')

View(RAD51B_Mutated_SV_per_Sample)

RAD51B_Mutated_SV_per_Sample <- data.table(RAD51B_Mutated_SV)[,.N, by=RAD51B_Mutated_SV$icgc_sample_id]
setnames(RAD51B_Mutated_SV_per_Sample, 'N', 'RAD51B Mutant SV Mutation Frequency')

median(RAD51B_Mutated_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`)
median(RAD51B_Not_Mutated_SV_per_Sample$`RAD51B WT SV Mutation Frequency`)



WT_boxplot <- ggplot(RAD51B_Not_Mutated_SV_per_Sample, aes(x=factor(0), y=RAD51B_Not_Mutated_SV_per_Sample$`RAD51B WT SV Mutation Frequency`))+geom_boxplot()+ylab('SVs per Sample')+xlab('RAD51B Not Mutated')+scale_y_log10(limits=c(1,100000))

MT_boxplot <- ggplot(RAD51B_Mutated_SV_per_Sample, aes(x=factor(0), y=RAD51B_Mutated_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`))+geom_boxplot()+ylab('SVs per Sample')+xlab('RAD51B Mutated')+scale_y_log10(limits=c(1,100000))


SVs_per_Sample_significance <- wilcox.test(RAD51B_Not_Mutated_SV_per_Sample$`RAD51B WT SV Mutation Frequency`,RAD51B_Mutated_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`, alternative = 'two.sided')$p.value

grid.arrange(WT_boxplot, MT_boxplot, ncol = 2, top = paste(('Wilcoxon, p<'), toString(SVs_per_Sample_significance)))



###Permutation Test

#NH1 - there is no difference in the mean SV per sample in RAD51B-ve samples compared to RAD51B-mutated samples
#AH1 - there is a difference in the mean SV per sample in RAD51B-ve samples compared to RAD51B-mutated samples

#NH2 - there is no difference in the median SV per sample in RAD51B-ve samples compared to RAD51B-mutated samples
#AH2 - there is a difference in the median SV per sample in RAD51B-ve samples compared to RAD51B-mutated samples


set.seed(1)

#Number of permutations
P <- 100000



Combined_Perm_Data <- setNames(data.frame(matrix(ncol=2, nrow=sum(nrow(RAD51B_Mutated_SV_per_Sample),nrow(RAD51B_Not_Mutated_SV_per_Sample)))), c('Mutation_Status','Sample_Mutation_Frequency'))

Combined_Perm_Data$Mutation_Status[1:nrow(RAD51B_Not_Mutated_SV_per_Sample)] <- 'WT'
Combined_Perm_Data$Mutation_Status[nrow(RAD51B_Not_Mutated_SV_per_Sample)+1:nrow(RAD51B_Mutated_SV_per_Sample)] <- 'Mutated'

for (i in 1:nrow(RAD51B_Not_Mutated_SV_per_Sample)) {
    Combined_Perm_Data$Sample_Mutation_Frequency[i] <- RAD51B_Not_Mutated_SV_per_Sample$`RAD51B WT SV Mutation Frequency`[i]
}


for (i in (nrow(RAD51B_Not_Mutated_SV_per_Sample)+1):(nrow(Combined_Perm_Data))) {
  Combined_Perm_Data$Sample_Mutation_Frequency[i] <- RAD51B_Mutated_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`[(i-402)]
}



#Number of observations to sample
n <- nrow(Combined_Perm_Data)

#Setting variable to resample from
sample_freq <- Combined_Perm_Data$Sample_Mutation_Frequency

#Calculate difference in sample means
mean_stat <- mean(Combined_Perm_Data$Sample_Mutation_Frequency[Combined_Perm_Data$Mutation_Status=='Mutated']) - mean(Combined_Perm_Data$Sample_Mutation_Frequency[Combined_Perm_Data$Mutation_Status=='WT'])
median_stat <- median(Combined_Perm_Data$Sample_Mutation_Frequency[Combined_Perm_Data$Mutation_Status=='Mutated']) - median(Combined_Perm_Data$Sample_Mutation_Frequency[Combined_Perm_Data$Mutation_Status=='WT'])


# initialize a matrix to store the permutation data
PermSamples <- data.frame(matrix(0, nrow=n, ncol=P))

PermSamples$Mutation_Status <- Combined_Perm_Data$Mutation_Status

View(PermSamples$Mutation_Status)

for(i in 2:P){
  PermSamples[,i] <- sample(sample_freq, size= n, replace=FALSE)
}

#PermSamples <- PermSamples[-c(2)]


#PermSamples <- PermSamples[,c(100001, 1:100000)]


Perm_mean_test_stat <- Perm_median_test_stat <- rep(0, P)

for (i in 2:P){
 
  Perm_mean_test_stat[i] <- mean(PermSamples[Combined_Perm_Data$Mutation_Status=="Mutated",i]) - 
                               mean(PermSamples[Combined_Perm_Data$Mutation_Status=="WT",i])

  Perm_median_test_stat[i] <- median(PermSamples[Combined_Perm_Data$Mutation_Status=="Mutated",i]) - 
                               median(PermSamples[Combined_Perm_Data$Mutation_Status=="WT",i])
}


Perm_mean_test_stat <- data.frame(Perm_mean_test_stat)
Perm_median_test_stat <- data.frame(Perm_median_test_stat)

#Perm_mean_test_stat <- Perm_mean_test_stat[-1,]
#Perm_median_test_stat <- Perm_median_test_stat[-1,]





#...calculate the p-value, for all P=100,000
sprintf(mean(Perm_mean_test_stat >= mean_stat), fmt = '%#.5g')
#Output 0.0000

# and, let's calculate the p-value for 
# option 2 of the test statistic (abs diff in medians)
sprintf(mean(Perm_median_test_stat >= median_stat), fmt = '%#.5g')
#Output 0.0000

#At 10%, 5% and 1% significance levels, we reject NH1 and NH2; there are significant differences in the SVs per sample when RAD51B is mutated compared to RAD51B-ve samples



## let's take a look at the 3 "Classic" hyp tests we could 
# consider (each of which comes with their own limitations...)

# Independent 2-sample t-test
# tests Ho: means are equal
t.test(Combined_Perm_Data$Sample_Mutation_Frequency ~ Combined_Perm_Data$Mutation_Status, paired=F, var.eq=F)  
#p-value = 7e-04


# let's look at the Wilcoxon aka Mann-Whitney U 
# tests Ho: medians are equal
wilcox.test(Combined_Perm_Data$Sample_Mutation_Frequency ~ Combined_Perm_Data$Mutation_Status, var.eq=F)  
#p-value<2e-16

# let's look at the Kolmogorov-Smirnov 2-sample test
# tests Ho: distributions are same
ks.test(Combined_Perm_Data$Sample_Mutation_Frequency[Combined_Perm_Data$Mutation_Status=='Mutated'],
        Combined_Perm_Data$Sample_Mutation_Frequency[Combined_Perm_Data$Mutation_Status=='WT'], paired = F)
#p-value = 1e-14


### Sampling distribution, and p-value for the Permutation approach, for mean values 

plot(density(Perm_mean_test_stat), 
     xlab=expression( group("|", bar(Yc) - bar(Ym), "|") ) , 
     main="Permutation Test Mean Distribution Plot", las=1, xlim = c(-80, 80))
abline(v=mean_stat, col="red", lty="dotted")
legend('topright', paste("p-value=",toString(sprintf(mean(Perm_mean_test_stat >= mean_stat), fmt = '%#.5g'))))



#Median values plot

plot(density(Perm_median_test_stat), 
     xlab=expression( group("|", bar(Yc) - bar(Ym), "|") ) , 
     main="Permutation Test Median Distribution Plot", las=1, xlim = c(-80, 80))
abline(v=median_stat, col="red", lty="dotted")
legend('topright', paste("p-value=",toString(sprintf(mean(Perm_median_test_stat >= median_stat), fmt = '%#.5g'))))


######
######


#Replicate Density Plot Analyses 

mean_SV_count <- mean(RAD51B_Mutated_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`)
median_SV_count <- median(RAD51B_Mutated_SV_per_Sample$`RAD51B Mutant SV Mutation Frequency`)


Perm_mean_SV_count <- Perm_median_SV_count <- rep(0, P)

for (i in 2:P){
  
  Perm_mean_SV_count[i] <- mean(PermSamples[Combined_Perm_Data$Mutation_Status=="Mutated",i]) 
   
  
  Perm_median_SV_count[i] <- median(PermSamples[Combined_Perm_Data$Mutation_Status=="Mutated",i])
  
}


Perm_mean_SV_count <- data.frame(Perm_mean_SV_count)
Perm_median_SV_count <- data.frame(Perm_median_SV_count)


Perm_mean_SV_count <- Perm_mean_SV_count[-1,]
Perm_median_SV_count <- Perm_median_SV_count[-1,]



plot(density(Perm_mean_SV_count), 
     xlab= 'Mean SV Count' , 
     main="Permutation Test Mean SV Count Distribution Plot", las=1, xlim =c(80,150))
abline(v=mean_SV_count, col="red", lty="dotted")
legend('topright', paste("p-value=",toString(sprintf(mean(Perm_mean_SV_count >= mean_SV_count), fmt = '%#.5g'))))


plot(density(Perm_median_SV_count), 
     xlab= 'Median SV Count' , 
     main="Permutation Test Median SV Count Distribution Plot", las=1, xlim=c(45,65))
abline(v=median_SV_count, col="red", lty="dotted")
legend('topright', paste("p-value=",toString(sprintf(mean(Perm_median_SV_count >= median_SV_count), fmt = '%#.5g'))))




###Repeating Analysis 

path <- https://docs.icgc.org/ap/v1/genes/ENSG00000182185/projects/
  



columns <- as.data.frame(colnames(RAD51B_Mutated_SV))
View(columns)


for (i in 1:ncol(RAD51B_Mutated_SV)){
  print(length(unique(RAD51B_Mutated_SV[,i])))
}


View(columns)
