
#############################################################
# Code to compute realtive gene expression of PR10 for all four Barley Plants in the Honours Project (Data used is from the qPCR on the #14th of November 2023)
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
library("pcr)

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/Ellie/")
#set working directory-Ellie CPU
#############################################################

##################################################################
#Import the data
#################################################################

#For tab delinated files only. row.name removes the 1,2,3, at the side when importing
GE_4<-read.delim("GE_Data_PR10UBI_14112023.txt", row.names=1)
GE_4

#make vector list for the genotypes
mapping_4 <- rep(c('Barke', 'SL17', 'SL52', 'GP'), each = 3)
mapping_4

##################################################################
#analyse the qPCR results
##################################################################

PCR_analysis_results_4_PR10 <- pcr_analyze(GE_4[-3], group_var = mapping_4, reference_gene = 'UBI', reference_group = 'SL52')
PCR_analysis_results_4_PR10

#results for this specifc qPCR
# group gene normalized calibrated relative_expression     error    lower    upper
#1 Barke PR10 -0.1654873 -1.0825672            2.117801 0.2886666 1.733756 2.586917
#2    GP PR10  0.6610330 -0.2560469            1.194202 0.4398076 0.880405 1.619844
#3  SL17 PR10 -0.1503245 -1.0674044            2.095660 0.5111932 1.470403 2.986794
#4  SL52 PR10  0.9170799  0.0000000            1.000000 0.1903091 0.876418 1.141008

##################################################################
#Plot Results of Gene Expression
##################################################################

#Genotype specific color which is color blind friendly

Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442")

PCR_Plot_4_PR10 <- pcr_analyze(GE_4[-3],group_var = mapping_4, reference_gene = 'UBI',reference_group = 'SL52', plot=TRUE) + labs(x = 'Barley Genotypes', y = 'Relative mRNA Expression') + ggtitle(label = 'PR10 Relative Gene Expression with SL52 Reference', subtitle='14th Nov. 2023') +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + geom_bar(stat="identity", fill=Genotype_Color, colour="black") + ylim(0,4.5)
PCR_Plot_4_PR10

##################################################################
#Determine Statistical Significance
##################################################################

#Linear regression used as there are more than 1 genotype
test_4_PR10 <-pcr_test(GE_4[-3], group_var = mapping_4, reference_gene = 'UBI', reference_group = 'SL52', test = 'lm')
test_4_PR10

#laid out in table
knitr::kable(test_4_PR10, caption = 'Figure 1: PR10 Linear Regression qPCR statistical Analysis')

#results of statistical analysis

#Table: Figure 1: PR10 Linear Regression qPCR statistical Analysis

#|gene |term           |   estimate|   p_value|      lower|      upper|
#|:----|:--------------|----------:|---------:|----------:|----------:|
#|PR10 |group_varBarke | -1.0825672| 0.0006684| -1.5473149| -0.6178195|
#|PR10 |group_varGP    | -0.2560469| 0.2396195| -0.7207946|  0.2087007|
#|PR10 |group_varSL17  | -1.0674044| 0.0007317| -1.5321521| -0.6026568|


##################################################################
#Plot Results of Gene Expression with Positive and Negative Standard Deviation
##################################################################

#Genotype specific color which is color blind friendly
Genotype_Color <- c("#0072B2", "#009E73","#56B4E9","#F0E442") 

#Relative Gene Expression Data taken from previously run 'pcr' Rstudio programme (line 51)
Relative_14_Data <-c(2.117801,1.194202, 2.095660, 1.000000)
Relative_14_Data

#standard deviation
sd_Relative_14 <- sd(Relative_14_Data)
sd_Relative_14

#Dataframe with all appropriate data 
PR10_14_Bargraph <- data.frame(Relative_GE = Relative_14_Data, Genotype = c("Barke", "GP", "SL17", "SL52"), sd= sd_Relative_14, yminvalue= Relative_14_Data-sd_Relative_14, ymaxvalue=Relative_14_Data+sd_Relative_14)
PR10_14_Bargraph

#plot
plot_PR10_14 <- ggplot(PR10_14_Bargraph,aes(x = Genotype, y = Relative_GE)) + geom_bar(stat = "identity", position = "dodge", colour="black") + labs(x = 'Barley Genotype', y = 'Relative Gene expression') + ggtitle(label = 'PR10 Relative Gene Expression with SL52 Reference', subtitle='14th Nov. 2023') +  geom_bar(stat="identity", fill=Genotype_Color, colour="black") +  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + ylim(-0.2,5.3) + geom_errorbar(aes(ymin=yminvalue, ymax=ymaxvalue, group=Genotype), width=0.4, colour='black', alpha=1, position = position_dodge(.9), linewidth=.5)

plot_PR10_14 #positive and negative standard deviation included