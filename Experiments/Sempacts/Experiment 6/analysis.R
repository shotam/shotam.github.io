##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Sempacts/'
exp = "Experiment 6/"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp6 = read.table(paste0(wd, exp, 'data.dat'), header=T)
Exp6 = recode.vars(Exp6, c("PrimeType", "Homogeneity"), contr.sum, scaleLevelDiffsTo1=T)

##### Create full model to compare against #####
Exp6.full<- glmer(TargBinary ~ PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev
                  + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                  , data = Exp6, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp6.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp6.NoPrimeType<- glmer(TargBinary ~ Homogeneity.lev + PrimeType.lev:Homogeneity.lev
                         + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                         , data = Exp6, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp6.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without Homogeneity #####
Exp6.NoHomogeneity<- glmer(TargBinary ~ PrimeType.lev + PrimeType.lev:Homogeneity.lev
                           + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                           , data = Exp6, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp6.NoHomogeneity@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x Homogeneity interaction #####

Exp6.NoInteraction<- glmer(TargBinary ~ PrimeType.lev + Homogeneity.lev 
                           + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                           , data = Exp6, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp6.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp6.full, Exp6.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_Homogeneity <- anova(Exp6.full, Exp6.NoHomogeneity); capture.output(MLR_Homogeneity, file = paste0(wd,exp, 'MLR_Homogeneity.txt'))
MLR_interaction <- anova(Exp6.full, Exp6.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))

# 
# > summary(Exp6.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBinary ~ PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev +      (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev |  
#                                                                                                 Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev |      Item)
# Data: Exp6
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 1889.1   2017.2   -920.5   1841.1     1512 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.3140 -0.7375 -0.2967  0.7875  3.0528 
# 
# Random effects:
#   Groups  Name                          Variance Std.Dev. Corr             
# Subject (Intercept)                   0.07222  0.2687                    
# PrimeType.lev                 0.46494  0.6819   -0.09            
# Homogeneity.lev               0.23426  0.4840    0.31  0.68      
# PrimeType.lev:Homogeneity.lev 0.09917  0.3149   -0.90 -0.20 -0.68
# Item    (Intercept)                   1.08522  1.0417                    
# PrimeType.lev                 0.40393  0.6356    0.33            
# Homogeneity.lev               0.20861  0.4567    0.90 -0.12      
# PrimeType.lev:Homogeneity.lev 0.63923  0.7995    0.34  1.00 -0.11
# Number of obs: 1536, groups:  Subject, 48; Item, 32
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    -0.3762     0.1994  -1.886  0.05923 .  
# PrimeType.lev                  -1.0970     0.1981  -5.537 3.08e-08 ***
#   Homogeneity.lev                -0.1510     0.1673  -0.903  0.36658    
# PrimeType.lev:Homogeneity.lev  -0.8226     0.2978  -2.763  0.00573 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. Hmgnt.
# PrimeTyp.lv 0.207               
# Homognty.lv 0.477  0.174        
# PrmTyp.l:H. 0.163  0.370  0.051 
