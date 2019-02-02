##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Sempacts/'
exp = "Experiment 2/"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp2 = read.table(paste0(wd, exp, 'data.dat'), header=T) 
Exp2 = recode.vars(Exp2, c("PrimeType","ExpCond","PrimeCategory"), contr.sum, scaleLevelDiffsTo1=T)
Exp2$PrimeType.lev = -1*Exp2$PrimeType.lev


##### Create full model to compare against #####
Exp2.full<- glmer(TargBinary ~ PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev
                       + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                       ,data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                
with(Exp2.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp2.NoPrimeType<- glmer(TargBinary ~ ExpCond.lev + PrimeType.lev:ExpCond.lev
                       + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                       , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                       
with(Exp2.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)
       
##### Create model without ExpCond #####
Exp2.NoExpCond<- glmer(TargBinary ~ PrimeType.lev + PrimeType.lev:ExpCond.lev
                       + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                       , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                       
with(Exp2.NoExpCond@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x ExpCond interaction #####

Exp2.NoInteraction<- glmer(TargBinary ~ PrimeType.lev + ExpCond.lev 
                       + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                       , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                  
with(Exp2.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp2.full, Exp2.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_ExpCond <- anova(Exp2.full, Exp2.NoExpCond); capture.output(MLR_ExpCond, file = paste0(wd,exp, 'MLR_ExpCond.txt'))
MLR_interaction <- anova(Exp2.full, Exp2.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))
write.csv(summary(Exp2.full)$coefficients, paste0(wd, exp, 'full_output.csv'))

# > summary(Exp2.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBinary ~ PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev +      (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |  
#                                                                                         Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |      Trial)
# Data: Exp2
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 3618.6   3765.8  -1785.3   3570.6     3381 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -8.8012 -0.6503  0.2746  0.6143  4.9685 
# 
# Random effects:
#   Groups  Name                      Variance Std.Dev. Corr             
# Subject (Intercept)               0.55910  0.7477                    
# PrimeType.lev             0.10008  0.3164   0.05             
# ExpCond.lev               0.05263  0.2294   0.13  1.00       
# PrimeType.lev:ExpCond.lev 0.16912  0.4112   0.82  0.61  0.67 
# Trial   (Intercept)               2.28850  1.5128                    
# PrimeType.lev             0.42908  0.6550   -0.72            
# ExpCond.lev               0.01416  0.1190    0.72 -0.03      
# PrimeType.lev:ExpCond.lev 0.06112  0.2472   -0.24 -0.50 -0.85
# Number of obs: 3405, groups:  Subject, 92; Trial, 48
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                0.90642    0.23864   3.798 0.000146 ***
#   PrimeType.lev             -0.82146    0.14328  -5.733 9.84e-09 ***
#   ExpCond.lev                0.13169    0.10212   1.290 0.197204    
# PrimeType.lev:ExpCond.lev -0.05871    0.20387  -0.288 0.773354    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. ExpCn.
# PrimeTyp.lv -0.480              
# ExpCond.lev  0.131  0.049       
# PrmTyp.:EC.  0.013 -0.017 -0.239
