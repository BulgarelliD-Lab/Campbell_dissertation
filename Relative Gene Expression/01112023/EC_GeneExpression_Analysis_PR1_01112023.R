
#############################################################
# Code to compute realtive gene expression of PR1 for all four Barley Plants in the Honours Project (Data used is from the qPCR on the #1st of November 2023)
#  Revision 12/23
#
#  Primary Contact Ellie Campbell: 2416433@dundee.ac.uk
#  Supervisor Dr. Davide Bulgarelli: d.bulgarelli@dundee.ac.uk 
#############################################################
#  Revision 12/23
#
#  Primary Contact Ellie Campbell: 2416433@dundee.ac.uk
#  Supervisor Dr. Davide Bulgarelli: d.bulgarelli@dundee.ac.uk 
#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#install once
install.packages("pcr")

#required packages
library("ggplot2")
library("PMCMR")
library("pcr")

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/Ellie/")
#set working directory-Ellie CPU
#############################################################

##################################################################
#Import the data
#################################################################

#For tab delinated files only. row.name removes the 1,2,3, at the side when importing
GE_2<-read.delim("GE_Data_PR1UBI_01112023.txt", row.names=1)
GE_2

#make vector list for the genotypes
mapping_2 <- rep(c('Barke', 'SL17', 'SL52', 'GP), each = 3)
mapping_2

##################################################################
#analyse the qPCR results
##################################################################

PCR_analysis_results_2_PR1 <- pcr_analyze(GE_2, group_var = mapping_2, reference_gene = 'UBI', reference_group = 'SL52', mode = 'seperate_tube')
PCR_analysis_results_2_PR1

#results for this specifc qPCR
#group gene normalized calibrated relative_expression     error       lower      upper
#1 Barke  PR1 -0.6912276  -6.841059          114.647313 0.1673207 102.0928481 128.745613
#2    GP  PR1  3.5311286  -2.618703            6.141975 0.3146578   4.9384037   7.638876
#3  SL17  PR1  1.3092689  -4.840562           28.651965 0.3662437  22.2281882  36.932164
#4  SL52  PR1  6.1498311   0.000000            1.000000 0.8160183   0.5680074   1.760540

##################################################################
#Plot Results of Gene Expression
##################################################################

#Genotype specific color which is color blind friendly

Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442")

PCR_Plot_2_PR1 <- pcr_analyze(GE_2,group_var = mapping_2, reference_gene = 'UBI',reference_group = 'SL52', plot=TRUE) + labs(x = 'Barley Genotypes', y = 'Relative mRNA expression') + ggtitle(label = 'PR1 Relative Gene Expression with SL52 Reference', subtitle='1st Nov. 2023') +   theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + geom_bar(stat="identity", fill=Genotype_Color, colour="black") + ylim(0,190)
PCR_Plot_2_PR1

##################################################################
#Determine Statistical Significance
##################################################################

#Linear regression used as there are more than 1 genotype 
test_2_PR1 <-pcr_test(GE_2, group_var = mapping_2, reference_gene = 'UBI', reference_group = 'SL52', test = 'lm')
test_2_PR1

#laid out in table
knitr::kable(test_2, caption = 'Figure 1: PR1 Linear Regression qPCR statistical Analysis')

#results of statistical analysis

#Table: Figure 1: PR1 Linear Regression qPCR statistical Analysis

#|gene |term           |  estimate|   p_value|     lower|     upper|
#|:----|:--------------|---------:|---------:|---------:|---------:|
#|PR1  |group_varBarke | -6.841059| 0.0000002| -7.805405| -5.876713|
#|PR1  |group_varGP    | -2.618703| 0.0002425| -3.583048| -1.654357|
#|PR1  |group_varSL17  | -4.840562| 0.0000028| -5.804908| -3.876216|


##################################################################
#Plot Results of Gene Expression with Positive and Negative Standard Deviation
##################################################################

#Genotype specific color which is color blind friendly
Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442") 

#Relative Gene Expression data taken from previously run 'pcr' Rstudio programme (line 51)
Relative_1_Data <-c(114.647313, 6.141975, 28.651965, 1.000000)
Relative_1_Data

#standard deviation
sd_Relative_1 <- sd(Relative_1_Data)
sd_Relative_1

#Dataframe with all appropriate data
PR10_1_Bargraph <- data.frame(Relative_GE = Relative_1_Data, Genotype = c("Barke", "GP", "SL17", "SL52"), sd= sd_Relative_1, yminvalue= Relative_1_Data-sd_Relative_1, ymaxvalue=Relative_1_Data+sd_Relative_1)
PR10_1_Bargraph

#plot
plot_PR1_1 <- ggplot(PR10_1_Bargraph,aes(x = Genotype, y = Relative_GE)) + geom_bar(stat = "identity", position = "dodge", colour="black") + labs(x = 'Barley Genotype', y = 'Relative Gene expression') + ggtitle(label = 'PR1 Relative Gene Expression with SL52 Reference', subtitle='1st Nov. 2023') +  geom_bar(stat="identity", fill=Genotype_Color, colour="black") +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + ylim(0,190) + geom_errorbar(aes(ymin=Relative_GE, ymax=ymaxvalue, group=Genotype), width=0.4, colour='black', alpha=1, position = position_dodge(.9), linewidth=.5)

plot_PR1_1 #only positive standard deviation included due large discrepancy between Barke and SL52