##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Sempacts/'
exp = "Experiment 5/"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp5 = read.table(paste0(wd, exp, 'data.dat'), header=T)
Exp5 = recode.vars(Exp5, c("PrimeType", "Homogeneity"), contr.sum, scaleLevelDiffsTo1=T)

##### Create full model to compare against #####
Exp5.full<- glmer(TargBinary ~ PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev
                  + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                  , data = Exp5, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp5.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp5.NoPrimeType<- glmer(TargBinary ~ Homogeneity.lev + PrimeType.lev:Homogeneity.lev
                         + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                         , data = Exp5, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp5.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without Homogeneity #####
Exp5.NoHomogeneity<- glmer(TargBinary ~ PrimeType.lev + PrimeType.lev:Homogeneity.lev
                       + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                       , data = Exp5, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp5.NoHomogeneity@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x Homogeneity interaction #####

Exp5.NoInteraction<- glmer(TargBinary ~ PrimeType.lev + Homogeneity.lev 
                           + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                           , data = Exp5, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp5.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp5.full, Exp5.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_Homogeneity <- anova(Exp5.full, Exp5.NoHomogeneity); capture.output(MLR_Homogeneity, file = paste0(wd,exp, 'MLR_Homogeneity.txt'))
MLR_interaction <- anova(Exp5.full, Exp5.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))

# > summary(Exp5.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBinary ~ PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev +      (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev |  
#                                                                                                 Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev |      Item)
# Data: Exp5
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 946.5   1063.2   -449.2    898.5      932 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -5.4799 -0.3731  0.2188  0.4329  2.8987 
# 
# Random effects:
#   Groups  Name                          Variance Std.Dev. Corr             
# Subject (Intercept)                   2.19036  1.4800                    
# PrimeType.lev                 0.17787  0.4217    0.50            
# Homogeneity.lev               1.49712  1.2236    0.16  0.94      
# PrimeType.lev:Homogeneity.lev 0.07668  0.2769   -0.10  0.81  0.97
# Item    (Intercept)                   2.36772  1.5387                    
# PrimeType.lev                 0.35950  0.5996   -0.10            
# Homogeneity.lev               1.23016  1.1091    0.61  0.72      
# PrimeType.lev:Homogeneity.lev 0.27339  0.5229    0.74  0.59  0.98
# Number of obs: 956, groups:  Subject, 48; Item, 45
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                     1.8142     0.3457   5.247 1.54e-07 ***
#   PrimeType.lev                  -0.3006     0.2911  -1.033    0.302    
# Homogeneity.lev                 0.5605     0.3650   1.536    0.125    
# PrimeType.lev:Homogeneity.lev  -0.2819     0.5379  -0.524    0.600    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. Hmgnt.
# PrimeTyp.lv 0.014               
# Homognty.lv 0.303  0.223        
# PrmTyp.l:H. 0.081  0.266  0.052
