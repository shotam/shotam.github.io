##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Pacts_OSF/'
exp = "Experiment3_"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp3 = read.csv(paste0(wd, exp, 'data.csv'), header=T)
Exp3 = recode.vars(Exp3, c("PrimeType", "Homogeneity"), contr.sum, scaleLevelDiffsTo1=T)
Exp3$PrimeType.lev = -1*Exp3$PrimeType.lev


### descriptive stats
table(Exp2$VerbClass)
Subj <- aggregate(data = Exp3, Response ~ Subject + VerbClass + Homogeneity + PrimeType, FUN = 'mean')
Overall <- aggregate(data = Subj, Response ~ PrimeType, FUN = 'mean')
Overall_byPrimeCat <- aggregate(data = Subj, Response ~ VerbClass + PrimeType, FUN = 'mean')
Overall_byHomogeneity <- aggregate(data = Subj, Response ~ Homogeneity + PrimeType, FUN = 'mean')

##### graphing #####

library(Rmisc)
Subj_noVerbClass <- aggregate(data = Exp3, Response ~ Subject+ Homogeneity + PrimeType, FUN = 'mean')
G <- summarySEwithin(data = Subj_noVerbClass, measurevar = 'Response', withinvars = c('PrimeType', 'Homogeneity'), idvar = 'Subject')
colnames(G) <- c("Prime_Structure", "Homogeneity", "n" , "mean", "sd", "se", "ci")

tabbedMeans <- tapply(G$mean, list(G$Prime_Structure, G$Homogeneity),
                      function(x) c(x = x))
tabbedSE <- tapply(G$se, list(G$Prime_Structure, G$Homogeneity),function(x) c(x = x))

pdf(paste0(wd,exp,'results.pdf'))

par(mar = c(5, 6, 4, 5) + 0.1)
barCenters <- barplot(height = tabbedMeans,
                      beside = TRUE, las = 1,
                      ylim = c(0, 1),
                      cex.lab = 1.2,
                      ylab = "Proportion of Target Response\nUsing the Preferred Syntactic Alternative",
                      xlab = "Event Homogeneity",
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

##### Create full model to compare against #####
Exp3.full<- glmer(Response ~ PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev
                  + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                  , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp3.NoPrimeType<- glmer(Response ~ Homogeneity.lev + PrimeType.lev:Homogeneity.lev
                         + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                         , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without Homogeneity #####
Exp3.NoHomogeneity<- glmer(Response ~ PrimeType.lev + PrimeType.lev:Homogeneity.lev
                       + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                       , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.NoHomogeneity@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x Homogeneity interaction #####

Exp3.NoInteraction<- glmer(Response ~ PrimeType.lev + Homogeneity.lev 
                           + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Subject) + (1 + PrimeType.lev + Homogeneity.lev + PrimeType.lev:Homogeneity.lev | Item)
                           , data = Exp3, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp3.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp3.full, Exp3.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_Homogeneity <- anova(Exp3.full, Exp3.NoHomogeneity); capture.output(MLR_Homogeneity, file = paste0(wd,exp, 'MLR_Homogeneity.txt'))
MLR_interaction <- anova(Exp3.full, Exp3.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))
