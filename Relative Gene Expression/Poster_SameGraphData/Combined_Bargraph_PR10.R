
#############################################################
# Code to compute plot relative gene expression results from PR10 on the same graph to be used in Honours Project Poster Presentation
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
library("pcr")

#############################################################
#set working directory-Davide CPU
setwd("/cluster/db/R_shared/Ellie/")
#set working directory-Ellie CPU
#############################################################

#############################################################
# Make a dataframe
#############################################################

#PR10 relative gene expression data taken from 'pcr' Rstudio programme: Code- EC_PR10_GeneExpression_Analysis_13112023 (line 51) & EC_GeneExpression_Analysis_PR10_141120231 (line 51)
#Plate 1= 13th Nov 2023, Plate 2= 14th Nov 2023

Relative_GE_Data_PR10 <-c(3.457057, 1.542498, 2.148518, 1.000000, 2.117801, 1.194202,  2.095660 ,1.000000 )
Relative_GE_Data_PR10

#STANDARD DEVIATION
Relative_GE_Data[0:4] #data Plate 1 13th November 2023 
sd_Relative_GE_1 <- sd(Relative_GE_Data_PR10[0:4])
sd_Relative_GE_1

Relative_GE_Data[4:8] #data Plate 2 14th November 2023
sd_Relative_GE_2 <- sd(Relative_GE_Data_PR10[5:8])
sd_Relative_GE_2

#Combine Standard Deviation
sd_Relative_GE <- c(rep(sd_Relative_GE_1, each= 4), rep(sd_Relative_GE_2, each= 4))
sd_Relative_GE

#Dataframe with all appropriate data
PR10_ALL_Bargraph <- data.frame(Relative_GE = Relative_GE_Data_PR10,  Plate_Type = rep(c('1', '2'), each = 4), Genotype = c("Barke", "GP", "SL17", "SL52"), sd= sd_Relative_GE, yminvalue= Relative_GE_Data_PR10-sd_Relative_GE, ymaxvalue=Relative_GE_Data_PR10+sd_Relative_GE)
PR10_ALL_Bargraph

# Relative_GE Plate_Type Genotype        sd yminvalue ymaxvalue
#1    3.457057          1    Barke 0.8251391 2.6319179  4.282196
#2    1.542498          1       GP 0.8251391 0.7173589  2.367637
#3    2.148518          1     SL17 0.8251391 1.3233789  2.973657
#4    1.000000          1     SL52 0.8251391 0.1748609  1.825139
#5    2.117801          2    Barke 0.8251391 1.2926619  2.942940
#6    1.194202          2       GP 0.8251391 0.3690629  2.019341
#7    2.095660          2     SL17 0.8251391 1.2705209  2.920799
#8    1.000000          2     SL52 0.8251391 0.1748609  1.825139

#############################################################
# Plot a Bargraph
#############################################################

#plot
plot_PR10_Both<- ggplot(PR10_ALL_Bargraph,aes(x = Plate_Type, y = Relative_GE, fill = Genotype )) + geom_bar(stat = "identity", position = "dodge", colour="black") + labs(x = 'Plate Replicate', y = 'Relative Gene expression') + ggtitle(label = 'PR10 Relative Gene Expression with SL52 Reference') +   theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 0.5)) + ylim(0,5.3) + scale_fill_manual(values=c("#0072B2", "#009E73","#56B4E9","#F0E442")) + geom_errorbar(aes(ymin=Relative_GE, ymax=ymaxvalue, group=Genotype), width=0.4, colour='black', alpha=1, position = position_dodge(.9), linewidth=.5)

plot_PR10_Both #Only positive standard error is graphed