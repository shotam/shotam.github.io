#### Tutorial for R ####

#### Data from semantic priming project: http://spp.montana.edu/
#### You can play around with full dataset, which contains insane amount of data.
#### In this tutorial, we will just focus on the first 48 subjects.
#df <- read.csv("/Users/shotamomma/Desktop/LING240_Rtutorial/data.csv")
#df <- subset(df, df$rel == 'rel' | df$rel == 'un')
#df_sub <- subset(df, df$subject <= 48)

#### load libraries ####
library(Rmisc)
library(lme4)
library(lmerTest)

#### load the data ####
df_sub <- read.csv("/Users/shotamomma/Desktop/LING240_Rtutorial/data_sub.csv")
df_correct <- subset(df_sub, df_sub$accuracy == 1) # you don't want to include the RT from incorrect trials, right?

#### make sure that R is interpreting each column correctly
df_correct$RT <- as.numeric(df_correct$RT)
df_correct$isi <- factor(df_correct$isi)
df_correct$rel <- factor(df_correct$rel)
df_correct$subject <- factor(df_correct$subject)
df_correct$target <- factor(df_correct$target)

#### coding the factors using deviation coding (very important when you interpret the model output)
#### R by default assigns 0 and 1. This does NOT allow you to get the estimate corresponding to ANOVA main effect!!!
contrasts(df_correct$isi) <- c(-0.5, 0.5)
contrasts(df_correct$rel) <- c(-0.5, 0.5)

#### trimming data

#### first trim the physically impossible RTs and RTs on the trials in which subjects zoned out.
df_correct$RT[df_correct$RT < 50] <- NA
df_correct$RT[df_correct$RT > 5000] <- NA

#### next trim the RTs that are more than 3 standard deviation from the mean for each participant.
df_correct$upper <- ave(df_correct$RT, df_correct$subject, FUN = function(x) mean(x, na.rm = T) + 3 * sd(x, na.rm = T))
df_correct$lower <- ave(df_correct$RT, df_correct$subject, FUN = function(x) mean(x, na.rm = T) - 3 * sd(x, na.rm = T))
df_correct$RT[df_correct$RT > df_correct$upper] <- NA
df_correct$RT[df_correct$RT < df_correct$lower] <- NA

#### let's take a look at F1 means and se ####
aggdf <- aggregate(df_correct$RT, list(df_correct$rel, df_correct$isi , df_correct$subject), mean, na.rm = T)
colnames(aggdf) <- c('relatedness', 'isi', 'subject', 'RT')
descriptiveStats <- summarySEwithin(aggdf, measurevar = 'RT', withinvars = c('relatedness', 'isi'), idvar = 'subject')

#relatedness  isi  N       RT       sd       se       ci
#1         rel   50 48 668.7350 47.75504 6.892846 13.86662
#2         rel 1050 48 683.1459 51.10437 7.376281 14.83916
#3          un   50 48 689.5818 54.56139 7.875258 15.84297
#4          un 1050 48 701.1615 50.97638 7.357806 14.80200

#### mixed effects model ####
m_min <- lmer(RT ~ rel*isi + (1|subject) + (1|target), data = df_correct) ### this is called 'intercept-only' model. Not acceptable in psycholinguistic
summary(m_min)

#Linear mixed model fit by REML 
#t-tests use  Satterthwaite approximations to degrees of freedom ['lmerMod']
#Formula: RT ~ rel * isi + (1 | subject) + (1 | target)
#Data: df_correct
#
#REML criterion at convergence: 505949.5
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-4.5552 -0.5730 -0.1673  0.3797  6.9354 
#
#Random effects:
#  Groups   Name        Variance Std.Dev.
#target   (Intercept)  2203     46.93  
#subject  (Intercept) 20329    142.58  
#Residual             34977    187.02  
#Number of obs: 37909, groups:  target, 1661; subject, 48

#Fixed effects:
#  Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)   686.085     20.634    47.000  33.250  < 2e-16 ***
#  rel1           19.708      1.923 36319.000  10.247  < 2e-16 ***
#  isi1           12.511      1.969 37847.000   6.355 2.11e-10 ***
#  rel1:isi1      -2.803      3.858 36624.000  -0.727    0.467    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Correlation of Fixed Effects:
#  (Intr) rel1   isi1  
#rel1       0.000              
#isi1       0.000 -0.002       
#rel1:isi1  0.000 -0.001  0.009

#### What is random effects? What is fixed effects? Don't ask me... ask someone reputable http://andrewgelman.com/2005/01/25/why_i_dont_use/

#### In general, you want to include every within-subject factor as by subject random slope, and every between-subject factor as by-item random slope
#### a model with maximal random effects structures. This can take a while to run, even with only 48 people's data
# m_max <- lmer(RT ~ rel*isi + (1 + rel*isi|subject) + (1 + rel*isi|target), data = df_correct) ### this is called a model with 'maximum random effect structure (Barr et al., 2013)'
#### a model with maximum random effects structure often doesn't converge...

#### when convergence failure occurs, remove the higest term random effects (or take a look at Brauer & Curtin (in press) table 17)
#m_1 <- lmer(RT ~ rel*isi + (1 + rel+isi|subject) + (1 + rel*isi|target), data = df_correct) ### this does not converge...
#m_2 <- lmer(RT ~ rel*isi + (1 + rel*isi|subject) + (1 + rel+isi|target), data = df_correct) ### Nope... 
m_3 <- lmer(RT ~ rel*isi + (1 + rel+isi|subject) + (1 + rel+isi|target), data = df_correct) ### YES!!!
# 
# > summary(m_3)
# Linear mixed model fit by REML 
# t-tests use  Satterthwaite approximations to degrees of freedom ['lmerMod']
# Formula: RT ~ rel * isi + (1 + rel + isi | subject) + (1 + rel + isi |      target)
# Data: df_correct
# 
# REML criterion at convergence: 504239
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.8264 -0.5708 -0.1670  0.3843  7.3250 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev. Corr       
# target   (Intercept)  2171.30  46.597             
# rel1          237.42  15.409  0.42       
# isi1           80.56   8.976  0.06  0.05 
# subject  (Intercept) 20360.34 142.690             
# rel1           23.73   4.871  -0.84      
# isi1         7043.35  83.925   0.20 -0.70
# Residual             33175.18 182.141             
# Number of obs: 37909, groups:  target, 1661; subject, 48
# 
# Fixed effects:
#   Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)   686.518     20.649    47.000  33.248   <2e-16 ***
#   rel1           19.856      2.037   190.000   9.750   <2e-16 ***
#   isi1           12.743     12.272    47.000   1.038    0.304    
# rel1:isi1      -3.300      3.765 34308.000  -0.876    0.381    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) rel1   isi1  
# rel1      -0.285              
# isi1       0.200 -0.238       
# rel1:isi1  0.000 -0.001  0.003

#### Once you buid the model, there are two ways to derive p-values ####
#### (1) maximum log liklihood ratio test
#### (2) lmerTest package (MCMC sample) - once you load lmerTest library, then summary() function will give you p-values
#### In the above example, we used lmerTest package
