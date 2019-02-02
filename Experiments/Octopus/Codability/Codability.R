library('Rmisc'); library('lme4');library('multcomp');library('ggplot2');library('lmerTest');library('ordinal');library(tidyr)
library('car')

df <- read.csv('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Codability/combined_data.csv')
df$participant <- interaction(df$participant, df$experiment)
df <- subset(df, df$participant != "211.Exp3")
df <- subset(df, df$participant != "15.Exp1")

df$participant <- factor(df$participant)
df$soa <- factor(df$soa)
df <- subset(df, df$Onset > 0)

df$Vcodability[df$Verb == 'bark'] <- 0.987
df$Vcodability[df$Verb == 'boil'] <- 2.223
df$Vcodability[df$Verb == 'bounce'] <- 1.645
df$Vcodability[df$Verb == 'cough'] <- 0.987
df$Vcodability[df$Verb == 'crawl'] <- 3.516
df$Vcodability[df$Verb == 'drown'] <- 1.807
df$Vcodability[df$Verb == 'fall'] <- 0.853
df$Vcodability[df$Verb == 'float'] <- 1.609
df$Vcodability[df$Verb == 'grow'] <- 1.567
df$Vcodability[df$Verb == 'melt'] <- 3.526
df$Vcodability[df$Verb == 'run'] <- 1.034
df$Vcodability[df$Verb == 'shake'] <- 2.612
df$Vcodability[df$Verb == 'shrink'] <- 1.311
df$Vcodability[df$Verb == 'sink'] <- 1.594
df$Vcodability[df$Verb == 'sleep'] <- 0.606
df$Vcodability[df$Verb == 'smile'] <- 2.516
df$Vcodability[df$Verb == 'sneeze'] <- 1.088
df$Vcodability[df$Verb == 'spin'] <- 2.874
df$Vcodability[df$Verb == 'swim'] <- 1.462
df$Vcodability[df$Verb == 'trip'] <- 3.001
df$Vcodability[df$Verb == 'walk'] <- 2.011
df$Vcodability[df$Verb == 'wink'] <- 1.513
df$Vcodability[df$Verb == 'yawn'] <- 0.771
df$Vcodability[df$Verb == 'yell'] <- 0.946


df$OnsetThe[df$OnsetThe > 5000 |df$OnsetThe < 300] <- NA ; sum(is.na(df$OnsetThe))/nrow(df)
df$Onset[df$Onset > 5000 |df$Onset < 300] <- NA; sum(is.na(df$Onset))/nrow(df)
df$The[df$The > 1500 |df$The < 30] <- NA; sum(is.na(df$The))/nrow(df)
df$dog[df$dog > 1500 |df$dog < 30 ] <- NA; sum(is.na(df$dog))/nrow(df)
df$above[df$above > 1500 |df$above < 30 ] <- NA; sum(is.na(df$above))/nrow(df)
df$the[df$the > 1500 |df$the < 30 ] <- NA;  sum(is.na(df$the))/nrow(df)
df$apple[df$apple > 1500 |df$apple < 30 ] <- NA;  sum(is.na(df$apple))/nrow(df)
df$is[df$is > 1500 |df$is < 30 ] <- NA;  sum(is.na(df$is))/nrow(df)
df$AI <- df$apple + df$is
df$AT <- df$above + df$the

##### log transform the data ####
df$OnsetThe <- log(df$OnsetThe)
df$Onset <- log(df$Onset)
df$The <- log(df$The)
df$dog <- log(df$dog)
df$above <- log(df$above)
df$the <- log(df$the)
df$apple <- log(df$apple)
df$is <- log(df$is)
df$AI <- log(df$AI)
df$AT <- log(df$AT)

##### sd trimming the data ####
sdvalue = 3
df$OnsetTheCritU <- ave(df$OnsetThe, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$OnsetTheCritL <- ave(df$OnsetThe, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$OnsetThe[df$OnsetTheCritU < df$OnsetThe | df$OnsetTheCritL > df$OnsetThe] <- NA

df$OnsetCritU <- ave(df$Onset, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$OnsetCritL <- ave(df$Onset, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$Onset[df$OnsetCritU < df$Onset | df$OnsetCritL > df$Onset] <- NA

df$TheCritU <- ave(df$The, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$TheCritL <- ave(df$The, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$The[df$TheCritU < df$The | df$TheCritL > df$The] <- NA

df$dogCritU <- ave(df$dog, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$dogCritL <- ave(df$dog, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$dog[df$dogCritU < df$dog | df$dogCritL > df$dog] <- NA

df$aboveCritU <- ave(df$above, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$aboveCritL <- ave(df$above, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$above[df$aboveCritU < df$above | df$aboveCritL > df$above] <- NA

df$theCritU <- ave(df$the, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$theCritL <- ave(df$the, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$the[df$theCritU < df$the | df$theCritL > df$the] <- NA

df$appleCritU <- ave(df$apple, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$appleCritL <- ave(df$apple, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$apple[df$appleCritU < df$apple | df$appleCritL > df$apple] <- NA

df$isCritU <- ave(df$is, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$isCritL <- ave(df$is, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$is[df$isCritU < df$is | df$isCritL > df$is] <- NA

df$AICritU <- ave(df$AI, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$AICritL <- ave(df$AI, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$AI[df$AICritU < df$AI | df$AICritL > df$AI] <- NA

df$ATCritU <- ave(df$AT, df$participant, FUN = function(x) mean(x, na.rm = T) + sdvalue * sd(x, na.rm = T))
df$ATCritL <- ave(df$AT, df$participant, FUN = function(x) mean(x, na.rm = T) - sdvalue * sd(x, na.rm = T))
df$AT[df$ATCritU < df$AT | df$ATCritL > df$AT] <- NA


df_unerg <- subset(df, df$VerbType == 'Unerg'); df_unerg$Verb <- factor(df_unerg$Verb)
tcodability <- tapply(df_unerg$Acodability, df_unerg$Verb, mean, na.rm = T)
t2 <- data.frame(tapply(df_unerg$OnsetThe, c(df_unerg$Verb), mean, na.rm = T))
t2$Vcodability <- tcodability
colnames(t2) <- c('rt', 'codability')


df_unacc <- subset(df, df$VerbType == 'Unacc'); df_unacc$Verb <- factor(df_unacc$Verb)
tcodability <- tapply(df_unacc$Acodability, df_unacc$Verb, mean, na.rm = T)
t2 <- data.frame(tapply(df_unacc$OnsetThe, c(df_unacc$Verb), mean, na.rm = T))
t2$Vcodability <- tcodability
colnames(t2) <- c('rt', 'codability')


library(rmcorr)
df_corr <- aggregate(df$OnsetThe, by = list(df$participant, df$Verb, df$Vcodability, df$VerbType), mean ,na.rm = T)
colnames(df_corr) <- c('participant', 'verb', 'Vcodability','VerbType', 'OnsetThe')
df_corr_unerg = subset(df_corr, df_corr$VerbType == 'Unacc')
codability.rmc <- rmcorr(participant = participant, measure1 = Vcodability, measure2 = OnsetThe, dataset = df_corr_unerg)
plot(codability.rmc, df_corr, overall = F, lty = 2, xlab = "x", ylab = "y")

df_corr <- aggregate(df$AI, by = list(df$participant, df$Verb, df$Vcodability, df$VerbType), mean ,na.rm = T)
colnames(df_corr) <- c('participant', 'verb', 'Vcodability','VerbType', 'OnsetThe')
df_corr_unerg = subset(df_corr, df_corr$VerbType == 'Unerg')
codability.rmc <- rmcorr(participant = participant, measure1 = Vcodability, measure2 = OnsetThe, dataset = df_corr_unerg)
plot(codability.rmc, df_corr, overall = F, lty = 2, xlab = "x", ylab = "y")


##### coding factors #####

df$VerbType = factor(df$VerbType); contrasts(df$VerbType) <- c(-0.5,0.5)  ### unacc = -0.5 ; unerg = 0.5
df$Relatedness = factor(df$Relatedness); contrasts(df$Relatedness) <- c(-0.5,0.5) ### related = -0.5 ; unrelated = 0.5
df$DistractorType = factor(df$DistractorType); contrasts(df$DistractorType) <- c(-0.5,0.5) ### 0.5 = Verb, -0.5 = noun
df$participant = factor(df$participant)
df$Transcription = factor(df$Transcription)
df$experiment = factor(df$experiment)

df$Vcodability <- scale(-df$Vcodability, center = T, scale = F)
df$Vagr <- scale(df$Vagr, center = T, scale = F)
df$Vfreq <- scale(log(df$Vfreq), center = T, scale = F)
df$Acodability <- scale(-df$Acodability, center = T, scale = F)
df$Afreq <- scale(log(df$Afreq), center = T, scale = F)
df$Scodability <- scale(-df$Scodability, center = T, scale = F)
df$Sfreq <- scale(log(df$Sfreq), center = T, scale = F)
df$unaccusativity <- scale(df$unaccusativity, center = T, scale = F)
df$VAoA <- scale(df$VAoA, center = T, scale = F)
df$Order <- scale(df$Order, center = T, scale = F)

df_region <- df %>% gather(Region, ProductionTime, Onset, The, dog, above, the, apple, is, OnsetThe, AI, AT)
df_region$Region <- factor(df_region$Region, ordered = T, levels = c("Onset", "The", "dog", "above", "the", "apple", "is","OnsetThe", "AI","AT"))

##### VerbDistractor #####

df_ROI <- subset(df_region, df_region$Region == 'OnsetThe'|df_region$Region == "AI")
df_ROI$Region <- factor(df_ROI$Region); contrasts(df_ROI$Region) <- c(-0.5,0.5)
df_ROI$participant <- factor(df_ROI$participant)

df_ROI_onset <- subset(df_ROI, df_ROI$Region == 'OnsetThe')
resid <- lmer(ProductionTime ~ unaccusativity + (1|participant),data = df_ROI_onset)
df_ROI_onset <- subset(df_ROI_onset, df_ROI_onset$ProductionTime > 0 & df_ROI_onset$unaccusativity != 'NA')
df_ROI_onset$resid <- residuals(resid)

codability_onset_model <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa*VerbType + Afreq + Afreq:VerbType + Acodability + Acodability:VerbType + Sfreq+ Sfreq:VerbType + Scodability:VerbType + Vfreq + Vfreq:VerbType + Vcodability + Vcodability:VerbType+ (1 + Vcodability*VerbType|participant) + (1|Verb), data = df_ROI_onset, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_onset_model_1 <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa*VerbType + Afreq + Afreq:VerbType + Acodability + Acodability:VerbType + Sfreq+ Sfreq:VerbType + Scodability:VerbType + Vfreq + Vfreq:VerbType + Vcodability+ (1 + Vcodability*VerbType|participant) + (1|Verb), data = df_ROI_onset, control=lmerControl(optCtrl=list(maxfun=100000)))

summary(codability_onset_model)

df_ROI_onset_unacc <- subset(df_ROI_onset, df_ROI_onset$VerbType == 'Unacc')
df_ROI_onset_unacc <- subset(df_ROI_onset_unacc, df_ROI_onset_unacc$ProductionTime > 0 & df_ROI_onset_unacc$unaccusativity != 'NA')
codability_onset_unacc_model <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq+ Vcodability+ + (1 + Vcodability|participant) + (1|Verb), data = df_ROI_onset_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_onset_unacc_model_1 <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq+ + (1 + Vcodability|participant) + (1|Verb), data = df_ROI_onset_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))

resid <- lmer(ProductionTime ~ unaccusativity + (1 + unaccusativity|participant) + (1|Verb),data = df_ROI_onset_unacc)
df_ROI_onset_unacc$resid <- residuals(resid)
codability_onset_unacc_model_resid <- lmer((resid) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq+ Vcodability+ + (1|participant) + (1|Verb), data = df_ROI_onset_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_onset_unacc_model_resid_1 <- lmer((resid) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq+ + (1|participant) + (1|Verb), data = df_ROI_onset_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))



r.squaredGLMM(codability_onset_unacc_model)

df_ROI_onset_unerg <- subset(df_region, df_region$Region == 'OnsetThe' & df_region$VerbType == 'Unerg')
codability_onset_unerg_model <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Scodability + Vcodability + Acodability  + (1 + Vcodability + Acodability|participant) + (1|Verb), data = df_ROI_onset_unerg, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_onset_unerg_model_1 <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Scodability + Acodability  + (1 + Vcodability + Acodability|participant) + (1|Verb), data = df_ROI_onset_unerg, control=lmerControl(optCtrl=list(maxfun=100000)))


df_ROI_AI <- subset(df_region, df_region$Region == 'AI')
codability_AI_model <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa*VerbType + Afreq + Afreq:VerbType + Acodability + Acodability:VerbType + Sfreq+ Sfreq:VerbType + Scodability + Scodability:VerbType + Vfreq + Vfreq:VerbType + Vcodability + Vcodability:VerbType+ (1 + Vcodability*VerbType|participant) + (1|Verb), data = df_ROI_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_AI_model_1 <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa*VerbType + Afreq + Afreq:VerbType + Acodability + Acodability:VerbType + Sfreq+ Sfreq:VerbType + Scodability + Scodability:VerbType + Vfreq + Vfreq:VerbType + Vcodability+ (1 + Vcodability*VerbType|participant) + (1|Verb), data = df_ROI_AI, control=lmerControl(optCtrl=list(maxfun=100000)))

df_ROI_AI_unacc <- subset(df_region, df_region$Region == 'AI' & df_region$VerbType == 'Unacc')
codability_AI_unacc_model <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq + Vcodability + (1 + Vcodability|participant) + (1|Verb), data = df_ROI_AI_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_AI_unacc_model_1 <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq + (1 + Vcodability|participant) + (1|Verb), data = df_ROI_AI_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))

df_ROI_AI_unerg <- subset(df_region, df_region$Region == 'AI' & df_region$VerbType == 'Unerg')
codability_AI_unerg_model <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq + Vcodability + (1 + Vcodability|participant) + (1|Verb), data = df_ROI_AI_unerg, control=lmerControl(optCtrl=list(maxfun=100000)))
codability_AI_unerg_model_1 <- lmer((ProductionTime) ~ Relatedness*DistractorType*soa + Afreq + Acodability + Sfreq+ Scodability + Vfreq + (1 + Vcodability|participant) + (1|Verb), data = df_ROI_AI_unerg, control=lmerControl(optCtrl=list(maxfun=100000)))

library(visreg)
pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Codability/EventCodability_ROI1.pdf')
visreg(codability_onset_model, "Vcodability", by = "VerbType", xlab = 'Event codability', ylab='Production Time (log)',
       points.par=list(col='black', cex=0.25), line = list(col = 'red'), gg = T) + theme_bw()
dev.off()

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Codability/EventCodability_ROI2.pdf')
visreg(codability_AI_model, "Vcodability", by = "VerbType", xlab = 'Event codability', ylab='Production Time (log)',
       points.par=list(col='black', cex=0.25), line = list(col = 'red'), gg = T) + theme_bw()
dev.off()


df_ROI <- subset(df_region, df_region$Region == 'OnsetThe'|df_region$Region == "AT")
#df_ROI$Region <- factor(df_ROI$Region, ordered = T, levels = c("OnsetThe", "AI")); contrasts(df_ROI$Region) <- c(-0.5,0.5)
df_ROI$Region <- factor(df_ROI$Region); contrasts(df_ROI$Region) <- c(-0.5,0.5)
df_ROI$participant <- factor(df_ROI$participant)

df_ROI_AT <- subset(df_region, df_region$Region == 'AT')
codability_AT_model <- lmer((ProductionTime) ~ Relatedness*VerbType*DistractorType*soa + Scodability + Scodability:VerbType + Sfreq + Sfreq:VerbType + Acodability + Acodability:VerbType +Afreq + Afreq:VerbType + Vcodability + Vcodability:VerbType + Vfreq + Vfreq:VerbType + (1 + Acodability + Acodability:VerbType |participant) + (1|Verb), data = df_ROI_AT, control=lmerControl(optCtrl=list(maxfun=100000)))


