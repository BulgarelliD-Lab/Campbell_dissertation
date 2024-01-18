
#############################################################
# Code to compute compare Ct values from two seperate runs of PR10 for all four Barley Plants in the Honours Project using 'pcr' package
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

#############################################################
# Import the Data sets
#############################################################

#For tab delinated files only. row.name removes the 1,2,3, at the side when importing
ALL_PR10 <- read.delim("Combined_PR10.txt", row.names=1)
ALL_PR10

#make vector list for the genotypes
Combined_Map<- rep(c('Barke', 'SL17', 'SL52', 'GP'), each=3)
Combined_Map

#############################################################
# Make a model matrix
#############################################################

group_PR10 <- relevel(factor(Combined_Map), refrence='SL52')
group_PR10

run <- factor(rep(c(1:2), each=12))
run

mm <- model.matrix(~group_PR10 + group_PR10:run, data = data.frame(group_PR10, run))
mm

#############################################################
# Run a Linear Regression to test for significant differences between runs
#############################################################

Combined_Test_PR10 <- pcr_test(ALL_PR10, reference_gene = 'UBI', model_matrix = mm, test = 'lm')
Combined_Test_PR10


knitr::kable(Combined_Test, caption = "Figure 19: Combining data from multiple qPCR runs")

##################################################################
#gene                            term   estimate      p_value     lower      upper
#1 PR10      model_matrixgroup_PR1Barke -1.7895444 0.0001199709 -2.541784 -1.0373044
#2 PR10         model_matrixgroup_PR1GP -0.6252683 0.0971459354 -1.377508  0.1269717
#3 PR10       model_matrixgroup_PR1SL17 -1.1033421 0.0067463633 -1.855582 -0.3511021
#4 PR10  model_matrixgroup_PR1SL52:run2 -1.2422784 0.0029579563 -1.994518 -0.4900384
#5 PR10 model_matrixgroup_PR1Barke:run2 -0.5353012 0.1509071225 -1.287541  0.2169388
#6 PR10    model_matrixgroup_PR1GP:run2 -0.8730571 0.0256342195 -1.625297 -0.1208171
#7 PR10  model_matrixgroup_PR1SL17:run2 -1.2063408 0.0036631216 -1.958581 -0.4541008
##################################################################



