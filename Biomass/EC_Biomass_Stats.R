
#############################################################
#  Code to compute calculations presented in the Honors_Project to analyse Biomass of Barley Plants
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
#required packages
library("ggplot2")
library("PMCMRplus")

#############################################################
#set working directory (If this doesn't work, do it via Session-> Select Working Directory -> Choose Directory)
setwd("/cluster/home/2416433/Honors_Project")

#############################################################

##################################################################
#Import the data
#################################################################

library(readxl)
GeneExpression_Biomass_R <- read_excel("GeneExpression_Biomass_R.xlsx")
View(GeneExpression_Biomass_R)
attach(GeneExpression_Biomass_R)

#inspect the file
GeneExpression_Biomass_R

###################################################################################################
#Step 1. let's define an experimental hypothesis
#Different genotypes have different Biomasses (g)

#Step 2 inspect independent variables and data distribution 
###################################################################################################

#identify independent variables
factor(Sample)
factor(Genotype)

#check dependent variables
sort(Weight)
hist(Weight)

#formal assessment of data distribution to identify the appropriate test 
#Shapiro test: https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test
shapiro.test(Weight)

###################################################################################################
#Shapiro-Wilk normality test

#data:  Weight
#W = 0.89077, p-value = 0.05729

#Interpret
# pvalue ≤ 0.05= population sample came from is NOT normally distriubted. Reject Null
# p-value> 0.05= population sample came from is normally distrubted. Accept Null
###################################################################################################

###################################################################################################
#Step 3 Prepare the datasets for and execute non parametric test
###################################################################################################

#Save each genotype separately 
Genotype_124_17 <- subset(GeneExpression_Biomass_R, Genotype == "124_17")
Genotype_124_52 <- subset(GeneExpression_Biomass_R, Genotype == "124_52")
Genotype_Barke <- subset(GeneExpression_Biomass_R, Genotype == "Barke")
Genotype_GP <- subset(GeneExpression_Biomass_R, Genotype == "GP")

#Inspect each dataset
Genotype_124_17 
Genotype_124_52 
Genotype_Barke 
Genotype_GP

#Basic Statistical Tests 

#MEAN
mean_124_17_Weight <- mean(Genotype_124_17$Weight) 
mean_124_52_Weight <- mean(Genotype_124_52$Weight)
mean_Barke_Weight <-mean(Genotype_Barke$Weight)
mean_GP_Weight <-mean(Genotype_GP$Weight)

#VARIANCE
variance_124_17_Weight <- var(Genotype_124_17$Weight) 
variance_124_52_Weight <- var(Genotype_124_52$Weight)
variance_Barke_Weight <-var(Genotype_Barke$Weight)
variance_GP_Weight <-var(Genotype_GP$Weight)

#Standard Deviation
sd_124_17_Weight <- sd(Genotype_124_17$Weight) 
sd_124_52_Weight <- sd(Genotype_124_52$Weight)
sd_Barke_Weight <- sd(Genotype_Barke$Weight)
sd_GP_Weight <- sd(Genotype_GP$Weight)

#Standard Error
se_124_17_Weight <- sd(Genotype_124_17$Weight)/sqrt(length(Genotype_124_17$Weight))
se_124_52_Weight <- sd(Genotype_124_52$Weight)/sqrt(length(Genotype_124_52$Weight))
se_Barke_Weight <- sd(Genotype_Barke$Weight)/sqrt(length(Genotype_Barke$Weight))
se_GP_Weight <- sd(Genotype_GP$Weight)/sqrt(length(Genotype_GP$Weight))

#########################################################################
#Effect of the Genotype on Biomass: If data is normally distributed (appropriate for this data set)
#########################################################################
aov(Weight ~ Genotype, data = GeneExpression_Biomass_R)
summary(aov(Weight ~ Genotype, data = GeneExpression_Biomass_R))

#########################################################################
#Call:
#  aov(formula = Weight ~ Genotype, data = GeneExpression_Biomass_R)

#Terms:
#                  Genotype  Residuals
#Sum of Squares  0.05377919 0.10588725
#Deg. of Freedom          3         12

#Residual standard error: 0.09393582
#Estimated effects may be unbalanced

#SUMMARY
#  Df  Sum Sq  Mean Sq F value Pr(>F)
#Genotype     3 0.05378 0.017926   2.032  0.163
#Residuals   12 0.10589 0.008824  

#########################################################################


#########################################################################
#Effect of the Genotype on Biomass: If data is not normally distributed (not appropriate for this data set)
#########################################################################

kruskal.test(Weight ~ Genotype, data = GeneExpression_Biomass_R)

#########################################################################
#RESULTS: THERE IS NO DIFFERENCE
#Kruskal-Wallis rank sum test
#data:  Weight by Genotype
#Kruskal-Wallis chi-squared = 5.5147, df = 3, p-value = 0.1378
#########################################################################

#post-hoc test: Pairwise comparison 

kwAllPairsDunnTest(GeneExpression_Biomass_R$Weight, as.factor(GeneExpression_Biomass_R$Genotype), p.adjust.method="BH")

#########################################################################
#Pairwise comparisons using Dunn's all-pairs test

#data: Weight and as.factor(Genotype)

#      124_17 124_52 Barke
#124_52 0.27   -      -    
#Barke  0.46   0.46   -    
#GP     0.16   0.46   0.27 

#P value adjustment method: BH
#alternative hypothesis: two.sided
########################################################################


###################################################################################################
#Step 4 Data visualisation
###########################################################################################

#Basic visualisation
p <- ggplot(GeneExpression_Biomass_R, aes(x=Genotype, y=Weight)) + geom_boxplot()
p

#something more elaborate...
EC_Genotype <- ordered(Genotype, levels=c("Barke", "GP", "124_17", "124_52"))
#So GP is golden Promise
levels(EC_Genotype) <- c("Barke", "Golden Promise", "SL17", "SL52")
#selected four colors
Genotype_Color_Biomass <- c("#0072B2", "#009E73","#56B4E9","#F0E442")

#Completed Plot with Legend Removed
Biomass_Plot<-ggplot(GeneExpression_Biomass_R, aes(x=EC_Genotype, y=Weight, color=EC_Genotype)) + ylab("Weight (g)") + xlab("Barley Genotype")+ geom_boxplot(lwd=1, outlier.shape = 17, outlier.size=2) + ggtitle('The Effect of Genotype on Biomass of Barley Plants') + geom_jitter(shape=16, position=position_jitter(0.2)) + scale_color_manual(values = Genotype_Color_Biomass) + labs(colour = "Genotype") + theme(plot.title = element_text(hjust = 0.5))+ theme(legend.position="none")

#Visualize Plot
Biomass_Plot

#Stats from geom_boxplot
layer_data(Biomass_Plot)

###########################################################################################
#Color can be linked to the appropriate genotype
#colour  ymin   lower middle   upper  ymax outliers notchupper notchlower x
#1 #0072B2 0.384 0.38400 0.4375 0.45575 0.473    0.261  0.4941825  0.3808175 1
#2 #009E73 0.164 0.25550 0.3075 0.35250 0.423           0.3841300  0.2308700 2
#3 #56B4E9 0.412 0.42400 0.4625 0.49775 0.500           0.5207625  0.4042375 3
#4 #F0E442 0.202 0.33025 0.3940 0.42500 0.455           0.4688525  0.3191475 4
###########################################################################################
