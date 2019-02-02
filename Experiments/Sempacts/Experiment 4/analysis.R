##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Sempacts/'
exp = "Experiment 4/"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp4 = read.table(paste0(wd, exp, 'data.dat'), header=T) 
Exp4 = recode.vars(Exp4, c("PrimeType", "CardTyp"), contr.sum, scaleLevelDiffsTo1=T)
Exp4$PrimeType.lev = -1*Exp4$PrimeType.lev

##### Create full model to compare against #####
Exp4.full<- glmer(TargBinary ~ PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev
                  + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Subject) + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Item)
                  , data = Exp4, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp4.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp4.NoPrimeType<- glmer(TargBinary ~ CardTyp.lev + PrimeType.lev:CardTyp.lev
                         + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Subject) + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Item)
                         , data = Exp4, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp4.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without CardTyp #####
Exp4.NoCardTyp<- glmer(TargBinary ~ PrimeType.lev + PrimeType.lev:CardTyp.lev
                       + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Subject) + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Item)
                       , data = Exp4, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp4.NoCardTyp@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x CardTyp interaction #####

Exp4.NoInteraction<- glmer(TargBinary ~ PrimeType.lev + CardTyp.lev 
                           + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Subject) + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev | Item)
                           , data = Exp4, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp4.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp4.full, Exp4.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_CardTyp <- anova(Exp4.full, Exp4.NoCardTyp); capture.output(MLR_CardTyp, file = paste0(wd,exp, 'MLR_CardTyp.txt'))
MLR_interaction <- anova(Exp4.full, Exp4.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))

# 
# > summary(Exp4.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBinary ~ PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev +      (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev |  
#                                                                                         Subject) + (1 + PrimeType.lev + CardTyp.lev + PrimeType.lev:CardTyp.lev |      Item)
# Data: Exp4
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 1757.5   1888.1   -854.8   1709.5     1683 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.3306 -0.6094  0.2409  0.5561  3.0218 
# 
# Random effects:
#   Groups  Name                      Variance Std.Dev. Corr             
# Subject (Intercept)               0.39188  0.6260                    
# PrimeType.lev             0.05644  0.2376    0.33            
# CardTyp.lev               0.03212  0.1792    0.77  0.86      
# PrimeType.lev:CardTyp.lev 0.20719  0.4552    0.91 -0.10  0.43
# Item    (Intercept)               2.97563  1.7250                    
# PrimeType.lev             0.44424  0.6665   -0.89            
# CardTyp.lev               0.13130  0.3624   -0.33  0.34      
# PrimeType.lev:CardTyp.lev 0.33681  0.5804    0.05  0.06 -0.91
# Number of obs: 1707, groups:  Subject, 48; Item, 48
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 1.3965     0.2844   4.910 9.09e-07 ***
#   PrimeType.lev              -1.1322     0.2085  -5.429 5.66e-08 ***
#   CardTyp.lev                -0.3199     0.1890  -1.693   0.0904 .  
# PrimeType.lev:CardTyp.lev   0.1130     0.3739   0.302   0.7624    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. CrdTy.
# PrimeTyp.lv -0.492              
# CardTyp.lev -0.109  0.202       
# PrmTyp.:CT.  0.112 -0.144 -0.466
