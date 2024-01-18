
#############################################################
# Code to compute realtive gene expression of PR10 for all four Barley Plants in the Honours Project (Data used is from the qPCR on the #13th of November 2023)
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
library("PMCMRplus")
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
GE_3_PR10<-read.delim("GE_Data_PR10UBI_13112023.txt", row.names=1)
GE_3_PR10

#make vector list for the genotypes
mapping_3 <- rep(c('Barke', 'SL17', 'SL52', 'GP'), each = 3)
mapping_3

##################################################################
#analyse the qPCR results
##################################################################

PCR_analysis_results_3_PR10 <- pcr_analyze(GE_3_PR10[-3], group_var = mapping_3, reference_gene = 'UBI', reference_group = 'SL52')
PCR_analysis_results_3_PR10

#results for this specifc qPCR
#group gene normalized calibrated relative_expression     error     lower    upper
#1 Barke PR10  0.3698139 -1.7895444            3.457057 0.2034663 3.0023208 3.980668
#2    GP PR10  1.5340900 -0.6252683            1.542498 0.5070935 1.0853609 2.192173
#3  SL17 PR10  1.0560163 -1.1033421            2.148518 0.5662250 1.4510699 3.181191
#4  SL52 PR10  2.1593583  0.0000000            1.000000 0.4544058 0.7298107 1.370218

##################################################################
#Plot Results of Gene Expression
##################################################################

#Genotype specific color which is color blind friendly

Genotype_Color <- c("#0072B2", "#56B4E9","#F0E442","#009E73")

PCR_Plot_3_PR10 <- pcr_analyze(GE_3_PR10[-3] ,group_var = mapping_3, reference_gene = 'UBI',reference_group = 'SL52', plot=TRUE) + labs(x = 'Barley Genotypes', y = 'Relative mRNA expression') + ggtitle(label = 'PR10 Relative Gene Expression with SL52 Reference', subtitle='13th Nov. 2023') +   theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + geom_bar(stat="identity", fill=Genotype_Color, colour="black") + ylim(0,4.5)
PCR_Plot_3_PR10

##################################################################
#Determine Statistical Significance
##################################################################

#Linear regression used as there are more than 1 genotype
test_3_PR10<-pcr_test(GE_3_PR10[-3], group_var = mapping_3, reference_gene = 'UBI', reference_group = 'SL52', test = 'lm')
test_3_PR10

#laid out in table
knitr::kable(test_3_PR10, caption = 'Figure 1: PR10 Linear Regression qPCR statistical Analysis')

#results of statistical analysis

#Table: Figure 1: PR10 Linear Regression qPCR statistical Analysis

#|gene |term           |   estimate|   p_value|     lower|      upper|
#|:----|:--------------|----------:|---------:|---------:|----------:|
#|PR10 |group_varBarke | -1.7895444| 0.0045841| -2.849338| -0.7297509|
#|PR10 |group_varGP    | -0.6252683| 0.2107564| -1.685062|  0.4345252|
#|PR10 |group_varSL17  | -1.1033421| 0.0431255| -2.163136| -0.0435486|


##################################################################
#Plot Results of Gene Expression with Positive and Negative Standard Deviation
##################################################################

#Genotype specific color which is color blind friendly
Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442") 

#Relative Gene Expression Data taken from previous ran 'pcr' Rstudio programme (line 51)
Relative_13_Data <-c(3.457057, 1.542498, 2.148518, 1.000000)
Relative_13_Data

#standard deviation
sd_Relative_13 <- sd(Relative_13_Data)
sd_Relative_13

#Data frame with appropriate data
PR10_13_Bargraph <- data.frame(Relative_GE = Relative_13_Data, Genotype = c("Barke", "GP", "SL17", "SL52"), sd= sd_Relative_13, yminvalue= Relative_13_Data-sd_Relative_13, ymaxvalue=Relative_13_Data+sd_Relative_13)
PR10_13_Bargraph

#plot
plot_PR10_13 <- ggplot(PR10_13_Bargraph,aes(x = Genotype, y = Relative_GE)) + geom_bar(stat = "identity", position = "dodge", colour="black") + labs(x = 'Barley Genotype', y = 'Relative Gene expression') + ggtitle(label = 'PR10 Relative Gene Expression with SL52 Reference', subtitle='13th Nov. 2023') +  geom_bar(stat="identity", fill=Genotype_Color, colour="black") +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + ylim(-0.2,5.3) + geom_errorbar(aes(ymin=yminvalue, ymax=ymaxvalue, group=Genotype), width=0.4, colour='black', alpha=1, position = position_dodge(.9), linewidth=.5)

plot_PR10_13 #positive and negative standard deviation included