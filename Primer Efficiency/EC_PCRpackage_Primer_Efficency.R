
#############################################################
#  Code to compute Primer Efficiency Honours Project using the 'pcr' package
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
PE_PR1UBIPR10<-read.delim("Primer_Analysis.txt", row.names=1)
PE_PR1UBIPR10

#make a vector of RNA amounts
amount <- rep(c(1, .3, .1, .03, .011, .0037, .0012), each = 3)
amount

##################################################################
#Standard Curve
#################################################################
#Determine R-Squared and Slope using Standard Curve Method(Slope required to calculate amplification efficiency)

standardcurve_UBIPR1PR10 <- pcr_assess(PE_PR1UBIPR10, amount = amount, method = 'standard_curve')
standardcurve_UBIPR1PR10

##################################################################
#gene intercept     slope r_squared
#1  UBI  20.69235 -3.448107 0.8137535
#2  PR1  23.56627 -3.101620 0.8970871
#3 PR10  22.78719 -3.040162 0.8864098

#Input slope into calculator For a graph where log (DNA copy#) is on the x-axis and Ct on the y-axis
#https://www.thermofisher.com/uk/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center #molecular-biology-resource-library/thermo-scientific-web-tools/qpcr-efficiency-calculator.html
##################################################################

##################################################################
#Plot with 'pcr' package
#################################################################

Standard_Plot_UBIPR1PR10 <- pcr_assess(PE_PR1UBIPR10, amount = amount, method = 'standard_curve', plot = TRUE) + labs(x = 'Quantity (Log(ng))', y = 'Ct value') + ggtitle(label = 'Primer Efficiency: Standard Curve', subtitle='Primers: UBI, PR1, PR10') +   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
Standard_Plot_UBIPR1PR10

##################################################################
#Plot with trendline and on individual plots
#################################################################
#Save slope and intercept. Below is specific for this dataset

#Interecept
UBI_Intercept<-20.69235
PR1_Intercept<-23.56627
PR10_Intercept<-22.78719

#Slope
UBI_Slope <- -3.448107
PR1_Slope <- -3.101620
PR10_Slope <--3.040162

#Import Specific Primer Data so can make individual plots with trendline

UBI_data<-read.delim("UBI_StandardCurve_Trendline.txt", row.names=1)
PR1_data<-read.delim("PR1_StandardCurve_Trendline.txt", row.names=1)
PR10_data<-read.delim("PR10_StandardCurve_Trendline.txt", row.names=1)


#Plot Primer Efficiency Curve with trendline

#Ubiquitin Primer Curve
UBI_primer_Plot <- plot(unlist(UBI_data)~log10(amount), xlab="Quantity (Log(ng))", ylab="Cycle Threshold (Ct)", main="Ubiquitin Primer Efficiency Curve") + abline(UBI_Intercept, UBI_Slope, col="red", lwd=2, lty=1)
UBI_primer_Plot

#PR1 Primer Curve
PR1_primer_Plot <- plot(unlist(PR1_data)~log10(amount), xlab="Quantity (Log(ng))", ylab="Cycle Threshold (Ct)", main="Pathogenesis Related 1 Primer Efficiency Curve") + abline(PR1_Intercept, PR1_Slope, col="red", lwd=2, lty=1)
PR1_primer_Plot

#PR10 Primer Curve
PR10_primer_Plot <- plot(unlist(PR10_data)~log10(amount), xlab="Quantity (Log(ng))", ylab="Cycle Threshold (Ct)", main="Pathogenesis Related 10 Primer Efficiency Curve") + abline(PR1_Intercept, PR1_Slope, col="red", lwd=2, lty=1)
PR10_primer_Plot