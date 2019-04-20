##### specify working directory #####
#wd = "/Volumes/LabFiles/Shota/Analysis/Sempacts/"
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Pacts_OSF/'
exp = "Experiment1/"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))


Exp1 = read.csv(paste0(wd, exp, 'data.csv'), header=T)
Exp1$PrimeCategory=factor(Exp1$PrimeCategory, levels = c("Dative","Transitive","Locative"))

Exp1 = recode.vars(Exp1, c("PrimeType", "PrimeCategory"), contr.sum, scaleLevelDiffsTo1=T)

### descriptive stats
Subj <- aggregate(data = Exp1, Response ~ Subject + PrimeCategory + PrimeType, FUN = 'mean')
Overall <- aggregate(data = Subj, Response ~ PrimeType, FUN = 'mean')
Overall_byPrimeCat <- aggregate(data = Subj, Response ~ PrimeCategory + PrimeType, FUN = 'mean')

##### graphing #####

library(Rmisc)

G <- summarySEwithin(data = Subj, measurevar = 'Response', withinvars = c('PrimeType', 'PrimeCategory'), idvar = 'Subject')
colnames(G) <- c("Prime_Structure", "Verb_Class", "n" , "mean", "sd", "se", "ci")
G$Verb_Class=factor(G$Verb_Class, levels = c("Dative","Locative","Transitive"))

tabbedMeans <- tapply(G$mean, list(G$Prime_Structure, G$Verb_Class),
                      function(x) c(x = x))
tabbedSE <- tapply(G$se, list(G$Prime_Structure, G$Verb_Class),function(x) c(x = x))

pdf(paste0(wd,exp,'results.pdf'))

par(mar = c(5, 6, 4, 5) + 0.1)
barCenters <- barplot(height = tabbedMeans,
                      beside = TRUE, las = 1,
                      ylim = c(0, 1),
                      cex.lab = 1.2,
                      ylab = "Proportion of Target Response\nUsing the Preferred Syntactic Alternative",
                      xlab = "Verb Class",
                      border = "black", axes = TRUE,
                      legend.text = TRUE,
                      args.legend = list(title = "Prime Struture", 
                                         xjust = 1, yjust = 1,
                                         x = "topright",
                                         cex = 1,
                                         bty ="n"))
segments(barCenters, tabbedMeans - tabbedSE, barCenters,
         tabbedMeans + tabbedSE, lwd = 1.5)

arrows(barCenters, tabbedMeans - tabbedSE, barCenters,
       tabbedMeans + tabbedSE, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

axis(1, at=c(0,10), labels = FALSE, tck=0)

dev.off()

##### Create full model to compare against: #####
Exp1.full<- glmer(Response ~ PrimeType.lev + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Item)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                       
#### including (1|PrimeCategory) or (1 + PrimeType.lev|PrimeCateogry) caused convergence failure.

summary(Exp1.full)
with(Exp1.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####

Exp1.NoPrimeType<- glmer(Response ~ PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Item) 
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
                       
summary(Exp1.NoPrimeType)
with(Exp1.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)  

##### Create model without PrimeCategory  #####
Exp1.NoPrimeCategory<- glmer(Response ~ PrimeType.lev + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Item)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))                       
                       
with(Exp1.NoPrimeCategory@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)                         
##### Create model without PrimeCategory x PrimeCategory interaction  #####
Exp1.NoInteraction<- glmer(Response ~ PrimeType.lev + PrimeCategory.lev1 + PrimeCategory.lev2 
                       + (1 + PrimeCategory.lev1 + PrimeCategory.lev2 + PrimeType.lev:PrimeCategory.lev1 + PrimeType.lev:PrimeCategory.lev2 | Subject) + (1 + PrimeType.lev | Item)
                       , data = Exp1, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp1.NoPrimeCategory@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)                         
##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp1.full, Exp1.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_PrimeCategory <- anova(Exp1.full, Exp1.NoPrimeCategory); capture.output(MLR_PrimeCategory, file = paste0(wd,exp, 'MLR_PrimeCategory.txt'))
MLR_interaction <- anova(Exp1.full, Exp1.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))
