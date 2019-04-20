##### specify working directory #####
wd = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Pacts_OSF/'
exp = "Experiment2_"

##### Pre-analysis #####
# Load packages + Kleiman re-corder
library(lme4)
source(paste0(wd, "Kleinman_recoder.R"))

# read data
Exp2 = read.csv(paste0(wd, exp, 'data.csv'), header=T) 
Exp2 = recode.vars(Exp2, c("PrimeType", "CardType"), contr.sum, scaleLevelDiffsTo1=T)
Exp2$PrimeType.lev = -1*Exp2$PrimeType.lev

### descriptive stats
table(Exp2$VerbClass)
Subj <- aggregate(data = Exp2, Response ~ Subject + VerbClass + CardType + PrimeType, FUN = 'mean')
Overall <- aggregate(data = Subj, Response ~ PrimeType, FUN = 'mean')
Overall_byCardType <- aggregate(data = Subj, Response ~ PrimeType + CardType, FUN = 'mean')
Overall_byPrimeCat <- aggregate(data = Subj, Response ~ VerbClass + CardType + PrimeType, FUN = 'mean')

##### graphing #####

library(Rmisc)

G <- summarySEwithin(data = Subj, measurevar = 'Response', withinvars = c('PrimeType', 'CardType'), idvar = 'Subject')
colnames(G) <- c("Prime_Structure", "CardType", "n" , "mean", "sd", "se", "ci")
G$CardType=factor(G$CardType, levels = c("Same","Different"))


tabbedMeans <- tapply(G$mean, list(G$Prime_Structure, G$CardType),
                      function(x) c(x = x))
tabbedSE <- tapply(G$se, list(G$Prime_Structure, G$CardType),function(x) c(x = x))

pdf(paste0(wd,exp,'results.pdf'))

par(mar = c(5, 6, 4, 5) + 0.1)
barCenters <- barplot(height = tabbedMeans,
                      beside = TRUE, las = 1,
                      ylim = c(0, 1),
                      cex.lab = 1.2,
                      ylab = "Proportion of Target Response\nUsing the Preferred Syntactic Alternative",
                      xlab = "Event Depiction",
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
Exp2.full<- glmer(Response ~ PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev
                  + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Subject) + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Item)
                  , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp2.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType #####
Exp2.NoPrimeType<- glmer(Response ~ CardType.lev + PrimeType.lev:CardType.lev
                         + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Subject) + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Item)
                         , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp2.NoPrimeType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without CardType #####
Exp2.NoCardType <- glmer(Response ~ PrimeType.lev + PrimeType.lev:CardType.lev
                       + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Subject) + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Item)
                       , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp2.NoCardType@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Create model without PrimeType x CardType interaction #####

Exp2.NoInteraction<- glmer(Response ~ PrimeType.lev + CardType.lev 
                           + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Subject) + (1 + PrimeType.lev + CardType.lev + PrimeType.lev:CardType.lev | Item)
                           , data = Exp2, family = binomial, verbose=T, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

with(Exp2.full@optinfo$derivs,max(abs(solve(Hessian,gradient)))<2e-3)

##### Conduct Maximum Likelihood Ratio tests #####
MLR_PrimeType <- anova(Exp2.full, Exp2.NoPrimeType); capture.output(MLR_PrimeType, file = paste0(wd,exp, 'MLR_PrimeType.txt'))
MLR_CardType <- anova(Exp2.full, Exp2.NoCardType); capture.output(MLR_CardTyp, file = paste0(wd,exp, 'MLR_CardType.txt'))
MLR_interaction <- anova(Exp2.full, Exp2.NoInteraction); capture.output(MLR_interaction, file = paste0(wd,exp, 'MLR_interaction.txt'))

