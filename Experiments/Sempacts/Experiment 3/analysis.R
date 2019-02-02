##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Sempacts/'
exp = "Experiment 3/"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp3 = read.table(paste0(wd, exp, 'data.dat'), header=T) 
Exp3 = recode.vars(Exp3, c("PrimeType", "ExpCond"), contr.sum, scaleLevelDiffsTo1=T)
Exp3$PrimeType.lev = -1*Exp3$PrimeType.lev

##### Create full model to compare against #####
Exp3.full<- glmer(TargBinary ~ PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev
                  + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                  , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp3.NoPrimeType<- glmer(TargBinary ~ ExpCond.lev + PrimeType.lev:ExpCond.lev
                         + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                         , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without ExpCond #####
Exp3.NoExpCond<- glmer(TargBinary ~ PrimeType.lev + PrimeType.lev:ExpCond.lev
                       + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                       , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.NoExpCond@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x ExpCond interaction #####

Exp3.NoInteraction<- glmer(TargBinary ~ PrimeType.lev + ExpCond.lev 
                           + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev | Trial)
                           , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp3.full, Exp3.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_ExpCond <- anova(Exp3.full, Exp3.NoExpCond); capture.output(MLR_ExpCond, file = paste0(wd,exp, 'MLR_ExpCond.txt'))
MLR_interaction <- anova(Exp3.full, Exp3.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))
write.csv(summary(Exp3.full)$coefficients, paste0(wd, exp, 'full_output.csv'))

# > summary(Exp3.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBinary ~ PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev +      (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |  
#                                                                                         Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |      Trial) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |      PrimeCategory)
# Data: Exp3
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 2173.3   2362.4  -1052.7   2105.3     1884 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -5.6078 -0.6963  0.2646  0.6517  2.1045 
# 
# Random effects:
#   Groups        Name                      Variance Std.Dev. Corr             
# Subject       (Intercept)               0.557076 0.74638                   
# PrimeType.lev             0.005544 0.07446  -0.31            
# ExpCond.lev               0.079591 0.28212  -0.64 -0.53      
# PrimeType.lev:ExpCond.lev 0.138613 0.37231  -0.64 -0.53  1.00
# Trial         (Intercept)               1.095023 1.04643                   
# PrimeType.lev             0.288481 0.53710   0.89            
# ExpCond.lev               0.007853 0.08862  -1.00 -0.88      
# PrimeType.lev:ExpCond.lev 0.059822 0.24459  -0.84 -0.50  0.85
# PrimeCategory (Intercept)               0.099620 0.31563                   
# PrimeType.lev             0.354650 0.59553   0.94            
# ExpCond.lev               0.027480 0.16577  -0.91 -0.71      
# PrimeType.lev:ExpCond.lev 0.020465 0.14306   0.45  0.72 -0.03
# Number of obs: 1918, groups:  Subject, 94; Trial, 24; PrimeCategory, 3
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                0.75346    0.29919   2.518   0.0118 *  
#   PrimeType.lev              1.53810    0.38641   3.981 6.88e-05 ***
#   ExpCond.lev               -0.04101    0.16829  -0.244   0.8075    
# PrimeType.lev:ExpCond.lev -0.57754    0.28688  -2.013   0.0441 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. ExpCn.
# PrimeTyp.lv  0.724              
# ExpCond.lev -0.446 -0.430       
# PrmTyp.:EC. -0.080  0.111  0.331


# 
# > summary(Exp3.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: TargBinary ~ PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev +      (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |  
#                                                                                         Subject) + (1 + PrimeType.lev + ExpCond.lev + PrimeType.lev:ExpCond.lev |      Trial)
# Data: Exp3
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 2159.9   2293.4  -1056.0   2111.9     1894 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -5.9908 -0.6904  0.2611  0.6604  2.0823 
# 
# Random effects:
#   Groups  Name                      Variance Std.Dev. Corr             
# Subject (Intercept)               0.558280 0.74718                   
# PrimeType.lev             0.005909 0.07687   0.36            
# ExpCond.lev               0.086126 0.29347  -0.65  0.48      
# PrimeType.lev:ExpCond.lev 0.146548 0.38282   0.62 -0.51 -1.00
# Trial   (Intercept)               1.193317 1.09239                   
# PrimeType.lev             0.669854 0.81845  -0.79            
# ExpCond.lev               0.023442 0.15311  -0.85  0.99      
# PrimeType.lev:ExpCond.lev 0.093042 0.30503   0.53  0.10 -0.01
# Number of obs: 1918, groups:  Subject, 94; Trial, 24
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                0.75475    0.24588   3.070  0.00214 ** 
#   PrimeType.lev             -1.55804    0.21702  -7.179 7.02e-13 ***
#   ExpCond.lev               -0.03394    0.14001  -0.242  0.80845    
# PrimeType.lev:ExpCond.lev  0.55801    0.27729   2.012  0.04418 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) PrmTy. ExpCn.
# PrimeTyp.lv -0.624              
# ExpCond.lev -0.255  0.261       
# PrmTyp.:EC.  0.175 -0.070 -0.397