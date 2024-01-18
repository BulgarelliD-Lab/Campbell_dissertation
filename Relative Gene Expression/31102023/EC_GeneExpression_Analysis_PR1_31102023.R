
#############################################################
# Code to compute realtive gene expression of PR1 for all four Barley Plants in the Honours Project (Data used is from the qPCR on the #31st of November 2023)
#  Revision 12/23
#
#  Primary Contact Ellie Campbell: 2416433@dundee.ac.uk
#  Supervisor Dr. Davide Bulgarelli: d.bulgarelli@dundee.ac.uk 
#############################################################
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

#import the datasets
#For tab delinated files only. row.name removes the 1,2,3, at the side when importing
GE_1<-read.delim("GE_Data_PR1UBI_31102023.txt", row.names=1)
GE_1

#make vector list for the genotypes
mapping_1 <- rep(c('Barke', 'SL17', 'SL52', 'GP'), each = 3)
mapping_1

##################################################################
#analyse the qPCR results
##################################################################

PCR_analysis_results_1_PR1 <- pcr_analyze(GE_1, group_var = mapping_1, reference_gene = 'UBI', reference_group = 'SL52')
PCR_analysis_results_1_PR1

#results for this specific qPCR
#group gene normalized calibrated relative_expression     error      lower      upper
#1 Barke  PR1 -0.1565329  -6.429866           86.214950 0.4735175 62.0925669 119.708654
#2    GP  PR1  4.5774892  -1.695844            3.239664 0.1834252  2.8528795   3.678887
#3  SL17  PR1  1.8668976  -4.406436           21.206515 0.2019760 18.4360747  24.393277
#4  SL52  PR1  6.2733332   0.000000            1.000000 0.1159127  0.9227983   1.083660

##################################################################
#Plot Results of Gene Expression
##################################################################

#Genotype specific color which is color blind friendly

Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442")

PCR_Plot_1_PR1 <- pcr_analyze(GE_1,group_var = mapping_1, reference_gene = 'UBI',reference_group = 'SL52', plot=TRUE) + labs(x = 'Barley Genotypes', y = 'Relative mRNA Expression') + ggtitle(label = 'PR1 Relative Gene Expression with SL52 Reference', subtitle='31st Oct. 2023') +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + geom_bar(stat="identity", fill=Genotype_Color, colour="black") + ylim(0,190)
PCR_Plot_1_PR1

##################################################################
#Determine Statistical Significance
##################################################################

#Linear regression used as there are more than 1 genotype 
test_1_PR1 <-pcr_test(GE_1, group_var = mapping_1, reference_gene = 'UBI', reference_group = 'SL52', test = 'lm')
test_1_PR1

#laid out in table
knitr::kable(test_1, caption = 'Figure 1: PR1 Linear Regression qPCR statistical Analysis')

#results of statistical analysis

#Table: Figure 1: PR1 Linear Regression qPCR statistical Analysis

#|gene |term           |  estimate|  p_value|     lower|     upper|
#|:----|:--------------|---------:|--------:|---------:|---------:|
#|PR1  |group_varBarke | -6.429866| 0.00e+00| -6.965894| -5.893839|
#|PR1  |group_varGP    | -1.695844| 8.42e-05| -2.231872| -1.159816|
#|PR1  |group_varSL17  | -4.406436| 1.00e-07| -4.942463| -3.870408|

##################################################################
#Plot Results of Gene Expression with Positive and Negative Standard Deviation
##################################################################

#Genotype specific color which is color blind friendly
Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442") 

#Relative Gene Expression data taken from previoulsy run 'pcr' Rstudio programme (line 51)
Relative_31_Data <-c(86.214950, 3.239664, 21.206515 ,1.000000)
Relative_31_Data

#standard deviation
sd_Relative_31 <- sd(Relative_31_Data)
sd_Relative_31

#Dataframe with all appropriate data
PR10_31_Bargraph <- data.frame(Relative_GE = Relative_31_Data, Genotype = c("Barke", "GP", "SL17", "SL52"), sd= sd_Relative_31, yminvalue= Relative_31_Data-sd_Relative_31, ymaxvalue=Relative_31_Data+sd_Relative_31)
PR10_31_Bargraph

#plot
plot_PR1_31 <- ggplot(PR10_1_Bargraph,aes(x = Genotype, y = Relative_GE)) + geom_bar(stat = "identity", position = "dodge", colour="black") + labs(x = 'Barley Genotype', y = 'Relative Gene expression') + ggtitle(label = 'PR1 Relative Gene Expression with SL52 Reference', subtitle='31st Oct. 2023') +  geom_bar(stat="identity", fill=Genotype_Color, colour="black") +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + ylim(0,190) + geom_errorbar(aes(ymin=Relative_GE, ymax=ymaxvalue, group=Genotype), width=0.4, colour='black', alpha=1, position = position_dodge(.9), linewidth=.5)

plot_PR1_31 #only positive standard deviation included due to large discrepancy between Barke and SL52