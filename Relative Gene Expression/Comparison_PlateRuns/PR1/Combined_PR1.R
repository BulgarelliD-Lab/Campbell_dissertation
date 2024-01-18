
#############################################################
# Code to compute compare Ct values from two seperate runs of PR1 for all four Barley Plants in the Honours Project using 'pcr' package
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
ALL_PR1 <- read.delim("Combined_PR1.txt", row.names=1)
ALL_PR1

#make vector list for the genotypes
Combined_Map<- rep(c('Barke', 'SL17', 'SL52', 'GP'), each=3)
Combined_Map

#############################################################
# Make a model matrix
#############################################################

group_PR1 <- relevel(factor(Combined_Map), refrence='SL52')
group_PR1

run <- factor(rep(c(1:2), each=12))
run

mm <- model.matrix(~group_PR1 + group_PR1:run, data = data.frame(group_PR1, run))
mm

#############################################################
# Run a Linear Regression to test for significant differences between runs
#############################################################

Combined_Test_PR1 <- pcr_test(ALL_PR1, reference_gene = 'UBI', model_matrix = mm, test = 'lm')
Combined_Test_PR1


knitr::kable(Combined_Test, caption = "Figure 19: Combining data from multiple qPCR runs")

##################################################################
#gene                            term   estimate      p_value      lower      upper
#1  PR1      model_matrixgroup_PR1Barke -6.4298662 2.097666e-12 -7.1470627 -5.7126696
#2  PR1         model_matrixgroup_PR1GP -1.6958440 1.275816e-04 -2.4130405 -0.9786475
#3  PR1       model_matrixgroup_PR1SL17 -4.4064356 6.222777e-10 -5.1236322 -3.6892391
#4  PR1  model_matrixgroup_PR1SL52:run2 -0.1235021 7.198545e-01 -0.8406986  0.5936944
#5  PR1 model_matrixgroup_PR1Barke:run2 -0.5346947 1.335651e-01 -1.2518912  0.1825018
#6  PR1    model_matrixgroup_PR1GP:run2 -1.0463607 6.983594e-03 -1.7635572 -0.3291642
#7  PR1  model_matrixgroup_PR1SL17:run2 -0.5576286 1.187951e-01 -1.2748251  0.1595679
##################################################################



