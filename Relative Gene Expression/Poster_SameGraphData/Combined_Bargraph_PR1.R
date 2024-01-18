
#############################################################
# Code to compute plot relative gene expression results from PR1 on the same graph to be used in Honours Project Poster Presentation
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

library("ggplot2")
library("PMCMRplus")

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/Ellie/")
#set working directory-Ellie CPU
#############################################################

#############################################################
# Make a dataframe
#############################################################

#PR10 relative gene expression data taken from 'pcr' Rstudio programme: Code- EC_GeneExpression_Analysis_PR1_31102023 (line 51) & EC_GeneExpression_Analysis_PR1_01112023 (line 51)
#Plate 1= 31st Oct 2023, Plate 2= 1st Nov 2023

Relative_GE_Data <-c(86.214950, 3.239664, 21.206515, 1.0000, 114.647313, 6.141975, 28.6519665, 1.0000)
Relative_GE_Data

#STANDARD DEVIATION
Relative_GE_Data[0:4] #data 31st October 2023
sd_Relative_GE_1 <- sd(Relative_GE_Data[0:4])
sd_Relative_GE_1

Relative_GE_Data[4:8] #data 1st November 2023
sd_Relative_GE_2 <- sd(Relative_GE_Data[5:8])
sd_Relative_GE_2

#Combined Standard Deviation
sd_Relative_GE <- c(rep(sd_Relative_GE_1, each= 4), rep(sd_Relative_GE_2, each= 4))
sd_Relative_GE

#Dataframe with all appropriate data
PR1_ALL_Bargraph <- data.frame(Relative_GE = Relative_GE_Data,  Plate_Type = rep(c('1', '2'), each = 4), Genotype = c("Barke", "GP", "SL17", "SL52"), sd= sd_Relative_GE, yminvalue= Relative_GE_Data-sd_Relative_GE, ymaxvalue=Relative_GE_Data+sd_Relative_GE)
PR1_ALL_Bargraph

#Relative_GE Plate_Type Genotype       sd yminvalue ymaxvalue
#1   86.214950          1    Barke 39.90479  46.31016 126.11974
#2    3.239664          1       GP 39.90479 -36.66513  43.14446
#3   21.206515          1     SL17 39.90479 -18.69828  61.11131
#4    1.000000          1     SL52 39.90479 -38.90479  40.90479
#5  114.647313          2    Barke 48.52261  66.12470 163.16992
#6    6.141975          2       GP 48.52261 -42.38063  54.66458
#7   28.651967          2     SL17 48.52261 -19.87064  77.17458
#8    1.000000          2     SL52 48.52261 -47.52261  49.52261

#############################################################
# Plot a Bargraph
#############################################################

#plot
plot_PR1_Both<- ggplot(PR1_ALL_Bargraph,aes(x = Plate_Type, y = Relative_GE, fill = Genotype )) + geom_bar(stat = "identity", position = "dodge", colour="black") + labs(x = 'Plate Replicate', y = 'Relative Gene expression') + ggtitle(label = 'PR1 Relative Gene Expression with SL52 Reference') +   theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + ylim(0,175) + scale_fill_manual(values=c("#0072B2", "#009E73","#56B4E9","#F0E442")) + geom_errorbar(aes(ymin=yminvalue, ymax=ymaxvalue, group=Genotype), width=0.4, colour='black', alpha=1, position = position_dodge(.9), linewidth=.5)

plot_PR1_Both #Only positive standard error is included due to large discrepency between SL52 and Barke