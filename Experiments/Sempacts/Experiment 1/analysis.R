##### specify working directory #####
#wd = "/Volumes/LabFiles/Shota/Analysis/Sempacts/"
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Sempacts/'
exp = "Experiment 1/"
##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp1 = read.table(paste0(wd,exp, 'data.dat'), header=T) 
Exp1 = recode.vars(Exp1, c("PrimeType", "PrimeCategory"), contr.sum, scaleLevelDiffsTo1=T)

##### Create full model to compare against: #####
Exp1.full<- glmer(TargBInary ~ PrimeType.lev + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Trial)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                       
#### including (1|PrimeCategory) or (1 + PrimeType.lev|PrimeCateogry) caused convergence failure.

summary(Exp1.full)
with(Exp1.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####

Exp1.NoPrimeType<- glmer(TargBInary ~ PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Trial) + (1+ PrimeType.lev|PrimeCategory)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                       
summary(Exp1.NoPrimeType)
with(Exp1.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)  

##### Create model without PrimeCategory  #####
Exp1.NoPrimeCategory<- glmer(TargBInary ~ PrimeType.lev + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Trial) + (1+ PrimeType.lev|PrimeCategory)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))                       
                       
with(Exp1.NoPrimeCategory@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)                         
##### Create model without PrimeCategory x PrimeCategory interaction  #####
Exp1.NoInteraction<- glmer(TargBInary ~ PrimeType.lev + PrimeCategory.lev1 + PrimeCategory.lev2 
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Trial) + (1+ PrimeType.lev|PrimeCategory)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp1.NoPrimeCategory@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)                         
##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp1.full, Exp1.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_PrimeCategory <- anova(Exp1.full, Exp1.NoPrimeCategory); capture.output(MLR_PrimeCategory, file = paste0(wd,exp, 'MLR_PrimeCategory.txt'))
MLR_interaction <- anova(Exp1.full, Exp1.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))
write.csv(summary(Exp1.full)$coefficients, paste0(wd, exp, 'full_output.csv'))
# 
# > summary(Exp1.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBInary ~ PrimeType.lev + PrimeCategory.lev1 + PrimeCategory.lev2 +      PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 +  
#   (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 +          PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev |      Trial)
# Data: Exp1
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 855.7    967.0   -403.8    807.7      742 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.5674 -0.4996  0.2146  0.5056  2.9464 
# 
# Random effects:
#   Groups  Name                             Variance Std.Dev. Corr                   
# Trial   (Intercept)                      2.4673   1.5708                          
# PrimeType.lev                    0.4166   0.6455   -0.44                  
# Subject (Intercept)                      0.8984   0.9479                          
# PrimeCategory.lev1               5.2128   2.2832   -0.69                  
# PrimeCategory.lev2               4.1006   2.0250    0.91 -0.55            
# PrimeCategory.lev1:PrimeType.lev 1.2084   1.0993    0.69 -0.99  0.50      
# PrimeCategory.lev2:PrimeType.lev 0.6311   0.7944   -0.49  0.96 -0.31 -0.97
# Number of obs: 766, groups:  Trial, 48; Subject, 22
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        0.9355     0.3377   2.770   0.0056 ** 
#   PrimeType.lev                     -1.1826     0.2782  -4.251 2.13e-05 ***
#   PrimeCategory.lev1                -0.4918     0.8722  -0.564   0.5729    
# PrimeCategory.lev2                -1.1892     0.8495  -1.400   0.1615    
# PrimeType.lev:PrimeCategory.lev1   0.4129     0.6976   0.592   0.5539    
# PrimeType.lev:PrimeCategory.lev2   0.6899     0.7082   0.974   0.3300    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. PrmC.1 PrmC.2 PT.:PC.1
# PrimeTyp.lv -0.278                              
# PrmCtgry.l1 -0.284  0.058                       
# PrmCtgry.l2  0.225  0.083 -0.464                
# PrmTy.:PC.1  0.205 -0.191 -0.379  0.130         
# PrmTy.:PC.2  0.011 -0.244  0.175 -0.263 -0.328  
# > 

