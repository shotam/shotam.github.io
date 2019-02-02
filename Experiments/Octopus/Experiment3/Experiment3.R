##### load data and packages #####

library('Rmisc'); library('lme4');library('multcomp');library('ggplot2');library('lmerTest');library('ordinal');library('tidyr');library('papaja')
df <- read.csv('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Codability/combined_data.csv')
df <- subset(df, df$experiment == "Exp4")
df$Picture<- interaction(df$Pic, df$Anoun)

sum(is.na(df$Onset))/nrow(df)
df <- subset(df, df$Onset > 0)

error <- data.frame(table(df$Relatedness, df$VerbType, df$soa))
error$percent <- error$Freq/(12*length(levels(as.factor(df$participant))))

df$OnsetThe <- df$Onset + df$The
df$OnsetTheRaw <- df$OnsetThe
df$dogRaw <- df$dog
df$aboveRaw <- df$above
df$theRaw <- df$the
df$appleRaw <- df$apple
df$isRaw <- df$is

##### exclude the extreme values #####
df$OnsetThe[df$OnsetThe > 5000 |df$OnsetThe < 30] <- NA ; sum(is.na(df$OnsetThe))/nrow(df)
df$Onset[df$Onset > 5000 |df$Onset < 30] <- NA; sum(is.na(df$Onset))/nrow(df)
df$The[df$The > 1500 |df$The < 30] <- NA; sum(is.na(df$The))/nrow(df)
df$dog[df$dog > 1500 |df$dog < 30 ] <- NA; sum(is.na(df$dog))/nrow(df)
df$above[df$above > 1500 |df$above < 30 ] <- NA; sum(is.na(df$above))/nrow(df)
df$the[df$the > 1500 |df$the < 30 ] <- NA;  sum(is.na(df$the))/nrow(df)
df$apple[df$apple > 1500 |df$apple < 30 ] <- NA;  sum(is.na(df$apple))/nrow(df)
df$is[df$is > 1500 |df$is < 30 ] <- NA;  sum(is.na(df$is))/nrow(df)

df$AI <- df$apple + df$is
df$AT <- df$above + df$the

df$AIRaw <- df$AI
df$ATRaw <- df$AT

df$pre <- (df$Onset + df$the + df$dog + df$above + df$the)/5
df$preN <- (df$Onset + df$the + df$dog)/3
df$post <- (df$dog + df$above + df$the + df$apple + df$is)/5

df$preRaw <- df$pre
df$preNRaw <- df$preN
df$postRaw <- df$post

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/distribution.pdf')
ggplot(subset(df, df$VerbType == "Unacc" & df$soa == "n150"), aes(OnsetThe, color=Relatedness)) +
  geom_histogram(position="identity", binwidth=50, aes(y=..density.., fill=Relatedness),  alpha=0.5) +
  geom_density()+
  theme_classic()
dev.off()

##### coding factors #####

df$VerbType = factor(df$VerbType); contrasts(df$VerbType) <- c(-0.5,0.5)  ### unacc = -0.5 ; unerg = 0.5
df$Relatedness = factor(df$Relatedness); contrasts(df$Relatedness) <- c(-0.5,0.5) ### related = -0.5 ; unrelated = 0.5
df$soa = factor(df$soa); contrasts(df$DistractorType) <- c(-0.5,0.5) ### 0.5 = Verb, -0.5 = noun
df$participant = factor(df$participant)
df$Transcription = factor(df$Transcription)

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
df$post <- log(df$post)
df$pre <- log(df$pre)
df$preN <- log(df$preN)

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


##### descriptive stats #####
df_region <- df %>% gather(Region, ProductionTime, Onset, The, dog, above, the, apple, is, OnsetThe ,AI, AT,pre,post)
df_region$Region <- factor(df_region$Region, ordered = T, levels = c("Onset", "The" , "dog", "above", "the", "apple", "is", "OnsetThe" , "AI","AT","pre","post"))

df_region_bySubj <- subset(df_region, df_region$Region != "OnsetThe" & df_region$Region != "AI" &df_region$Region != "AT" & df_region$Region != "pre" & df_region$Region != "post") %>% dplyr::group_by(participant,VerbType, Relatedness, soa, Region) %>%
  dplyr::summarise(ProductionTime = mean(exp(ProductionTime), na.rm = TRUE))

descriptives <- df_region_bySubj %>% dplyr::group_by(soa, VerbType, Relatedness, Region) %>%
  dplyr::summarise(
    Mean = as.integer(mean((ProductionTime), na.rm = TRUE)),
  )

descriptives <- spread(descriptives, Region, Mean)
apa_table(descriptives)


##### creating graph: region-by-region difference #####
dfNeg_graph <- subset(df, df$soa == 'n150')
graph1 <- as.data.frame(tapply(exp(dfNeg_graph$Onset),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph1$UnaccDiff <- graph1$Related.Unacc - graph1$Unrelated.Unacc
graph1$UnergDiff <- graph1$Related.Unerg - graph1$Unrelated.Unerg
graph1$region = '1'
graph2 <- as.data.frame(tapply(exp(dfNeg_graph$The),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph2$UnaccDiff <- graph2$Related.Unacc - graph2$Unrelated.Unacc
graph2$UnergDiff <- graph2$Related.Unerg - graph2$Unrelated.Unerg
graph2$region = '2'
graph3 <- as.data.frame(tapply(exp(dfNeg_graph$dog),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph3$UnaccDiff <- graph3$Related.Unacc - graph3$Unrelated.Unacc
graph3$UnergDiff <- graph3$Related.Unerg - graph3$Unrelated.Unerg
graph3$region = '3'
graph4 <- as.data.frame(tapply(exp(dfNeg_graph$above),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph4$UnaccDiff <- graph4$Related.Unacc - graph4$Unrelated.Unacc
graph4$UnergDiff <- graph4$Related.Unerg - graph4$Unrelated.Unerg
graph4$region = '4'
graph5 <- as.data.frame(tapply(exp(dfNeg_graph$the),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph5$UnaccDiff <- graph5$Related.Unacc - graph5$Unrelated.Unacc
graph5$UnergDiff <- graph5$Related.Unerg - graph5$Unrelated.Unerg
graph5$region = '5'
graph6 <- as.data.frame(tapply(exp(dfNeg_graph$apple),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph6$UnaccDiff <- graph6$Related.Unacc - graph6$Unrelated.Unacc
graph6$UnergDiff <- graph6$Related.Unerg - graph6$Unrelated.Unerg
graph6$region = '6'
graph7 <- as.data.frame(tapply(exp(dfNeg_graph$is),list(dfNeg_graph$participant,dfNeg_graph$Relatedness,dfNeg_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph7$UnaccDiff <- graph7$Related.Unacc - graph7$Unrelated.Unacc
graph7$UnergDiff <- graph7$Related.Unerg - graph7$Unrelated.Unerg
graph7$region = '7'
graph <- rbind(graph1,graph2,graph3,graph4,graph5,graph6, graph7)
graph$region <- factor(graph$region, ordered = T, levels = c('1','2','3','4','5','6','7'))

UnaccMean <- aggregate(graph$UnaccDiff, by= list(graph$region), FUN = function(x) mean(x, na.rm = T))
colnames(UnaccMean) <- c('region', 'mean')
UnaccSE <- aggregate(graph$UnaccDiff, by= list(graph$region), FUN = function(x) sd(x, na.rm = T)/sqrt(length(graph[,1])/6))
colnames(UnaccSE) <- c('region', 'se')
Unacc <- UnaccMean
Unacc$se <- UnaccSE$se
Unacc$VerbType <- 'Unacc'

UnergMean <- aggregate(graph$UnergDiff, by= list(graph$region), FUN = function(x) mean(x, na.rm = T))
colnames(UnergMean) <- c('region', 'mean')
UnergSE <- aggregate(graph$UnergDiff, by= list(graph$region), FUN = function(x) sd(x, na.rm = T)/sqrt(length(graph[,1])/6))
colnames(UnergSE) <- c('region', 'se')
Unerg <- UnergMean
Unerg$se <- UnergSE$se
Unerg$VerbType <- 'Unerg'

Graph <- rbind(Unacc, Unerg)

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/Negdifference.pdf')
ggplot(Graph, aes(x= region, y=mean, group=VerbType, color=VerbType)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05)) +
  
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7),labels=c("1" = "Onset", "2" = "the", "3" = "octopus",
                                                      "4" = "below", "5" = "the",
                                                      "6" = "spoon", "7" = "is")) +
  labs(title="Verb interference \n", x="\n Region", y = "Mean interference (ms) \n") +
  geom_hline(aes(yintercept=0)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.position = c(0.8, 0.8),
        axis.text.x  = element_text(angle=45, hjust = 1, size=14,face=c("bold","bold","italic","italic","italic","bold","bold"))) +
  coord_cartesian(ylim = c(-50, 100))
dev.off()


dfSim_graph <- subset(df, df$soa == 'pn0')
graph1 <- as.data.frame(tapply(exp(dfSim_graph$Onset),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph1$UnaccDiff <- graph1$Related.Unacc - graph1$Unrelated.Unacc
graph1$UnergDiff <- graph1$Related.Unerg - graph1$Unrelated.Unerg
graph1$region = '1'
graph2 <- as.data.frame(tapply(exp(dfSim_graph$The),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph2$UnaccDiff <- graph2$Related.Unacc - graph2$Unrelated.Unacc
graph2$UnergDiff <- graph2$Related.Unerg - graph2$Unrelated.Unerg
graph2$region = '2'
graph3 <- as.data.frame(tapply(exp(dfSim_graph$dog),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph3$UnaccDiff <- graph3$Related.Unacc - graph3$Unrelated.Unacc
graph3$UnergDiff <- graph3$Related.Unerg - graph3$Unrelated.Unerg
graph3$region = '3'
graph4 <- as.data.frame(tapply(exp(dfSim_graph$above),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph4$UnaccDiff <- graph4$Related.Unacc - graph4$Unrelated.Unacc
graph4$UnergDiff <- graph4$Related.Unerg - graph4$Unrelated.Unerg
graph4$region = '4'
graph5 <- as.data.frame(tapply(exp(dfSim_graph$the),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph5$UnaccDiff <- graph5$Related.Unacc - graph5$Unrelated.Unacc
graph5$UnergDiff <- graph5$Related.Unerg - graph5$Unrelated.Unerg
graph5$region = '5'
graph6 <- as.data.frame(tapply(exp(dfSim_graph$apple),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph6$UnaccDiff <- graph6$Related.Unacc - graph6$Unrelated.Unacc
graph6$UnergDiff <- graph6$Related.Unerg - graph6$Unrelated.Unerg
graph6$region = '6'
graph7 <- as.data.frame(tapply(exp(dfSim_graph$is),list(dfSim_graph$participant,dfSim_graph$Relatedness,dfSim_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph7$UnaccDiff <- graph7$Related.Unacc - graph7$Unrelated.Unacc
graph7$UnergDiff <- graph7$Related.Unerg - graph7$Unrelated.Unerg
graph7$region = '7'
graph <- rbind(graph1,graph2,graph3,graph4,graph5,graph6, graph7)
graph$region <- factor(graph$region, ordered = T, levels = c('1','2','3','4','5','6','7'))

UnaccMean <- aggregate(graph$UnaccDiff, by= list(graph$region), FUN = function(x) mean(x, na.rm = T))
colnames(UnaccMean) <- c('region', 'mean')
UnaccSE <- aggregate(graph$UnaccDiff, by= list(graph$region), FUN = function(x) sd(x, na.rm = T)/sqrt(length(graph[,1])/6))
colnames(UnaccSE) <- c('region', 'se')
Unacc <- UnaccMean
Unacc$se <- UnaccSE$se
Unacc$VerbType <- 'Unacc'

UnergMean <- aggregate(graph$UnergDiff, by= list(graph$region), FUN = function(x) mean(x, na.rm = T))
colnames(UnergMean) <- c('region', 'mean')
UnergSE <- aggregate(graph$UnergDiff, by= list(graph$region), FUN = function(x) sd(x, na.rm = T)/sqrt(length(graph[,1])/6))
colnames(UnergSE) <- c('region', 'se')
Unerg <- UnergMean
Unerg$se <- UnergSE$se
Unerg$VerbType <- 'Unerg'

Graph <- rbind(Unacc, Unerg)

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/Simdifference.pdf')
ggplot(Graph, aes(x= region, y=mean, group=VerbType, color=VerbType)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05)) +
  
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7),labels=c("1" = "Onset", "2" = "the", "3" = "octopus",
                                                      "4" = "below", "5" = "the",
                                                      "6" = "spoon", "7" = "is")) +
  labs(title="Noun interference \n", x="\n Region", y = "Mean interference (ms) \n") +
  geom_hline(aes(yintercept=0)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), legend.position = c(0.8, 0.8),
        axis.text.x  = element_text(angle=45, hjust = 1, size=14,face=c("bold","bold","italic","italic","italic","bold","bold"))) +
  coord_cartesian(ylim = c(-50, 100))
dev.off()

##### creating graphs: ROI analysis #####
library(Hmisc)
pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/NegOnset.pdf')
apa_lineplot(
  data = dfNeg_graph,
  id = "participant",
  dv = "OnsetThe",
  factors = c("VerbType", "Relatedness"),
  dispersion = wsci,
  las = 1,
  args_points = list(cex = 1),
  args_arrows = list(length = .075),
  args_y_axis = NULL,
  jitter = 0,
  ylim = c(7.1,7.3),
  ylab = 'Production Time (log)'
)
dev.off()

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/NegAI.pdf')
apa_lineplot(
  dfNeg_graph,
  id = "participant",
  dv = "AI",
  factors = c("VerbType", "Relatedness"),
  dispersion = wsci,
  las = 1,
  args_points = list(cex = 1),
  args_arrows = list(length = .075),
  args_y_axis = NULL,
  jitter = 0,
  ylim = c(6.55,6.75),
  ylab = 'Production Time (log)'
)
dev.off()

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/SimOnset.pdf')
apa_lineplot(
  dfSim_graph,
  id = "participant",
  dv = "OnsetThe",
  factors = c("VerbType", "Relatedness"),
  dispersion = wsci,
  las = 1,
  args_points = list(cex = 1),
  args_arrows = list(length = .075),
  args_y_axis = NULL,
  jitter = 0,
  ylim = c(7.10,7.30),
  ylab = 'Production Time (log)'
)
dev.off()

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/SimAT.pdf')
apa_lineplot(
  dfSim_graph,
  id = "participant",
  dv = "AT",
  factors = c("VerbType", "Relatedness"),
  dispersion = wsci,
  las = 1,
  args_points = list(cex = 1),
  args_arrows = list(length = .075),
  args_y_axis = NULL,
  jitter = 0,
  ylim = c(6.2,6.4),
  ylab = 'Production Time (log)'
)
dev.off()

##### Negative SOA #####

df_ROI <- subset(df_region, df_region$Region == 'OnsetThe'|df_region$Region == "post")
df_ROI$Region <- factor(df_ROI$Region, ordered = T, levels = c("OnsetThe", "post")); contrasts(df_ROI$Region) <- c(-0.5,0.5)

dfNeg <- subset(df_ROI, df_ROI$soa == 'n150')
dfSim <- subset(df_ROI, df_ROI$soa == 'pn0')

##### First ROI

dfNeg_OnsetThe <- subset(dfNeg, dfNeg$Region == "OnsetThe")
Neg_OnsetThe_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfNeg_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
Neg_OnsetThe_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfNeg_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
Neg_OnsetThe_inverse <- glmer((OnsetTheRaw) ~ Relatedness*VerbType+ (1 + Relatedness*VerbType - Relatedness - VerbType||participant) + (1+Relatedness||Verb), data = subset(df, df$soa == 'n150'), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) 

saveRDS(Neg_OnsetThe_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_log.rds")
saveRDS(Neg_OnsetThe_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_raw.rds")
saveRDS(Neg_OnsetThe_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_inverse.rds")


Neg_OnsetThe_log_coefs <- data.frame(coef(summary(Neg_OnsetThe_log)))
Neg_OnsetThe_log_coefs$p.z <- 2 * (1 - pnorm(abs(Neg_OnsetThe_log_coefs$t.value)))
Neg_OnsetThe_log_coefs
Neg_OnsetThe_raw_coefs <- data.frame(coef(summary(Neg_OnsetThe_raw)))
Neg_OnsetThe_raw_coefs$p.z <- 2 * (1 - pnorm(abs(Neg_OnsetThe_raw_coefs$t.value)))
Neg_OnsetThe_raw_coefs


dfSim_OnsetThe <- subset(dfSim, dfSim$Region == "OnsetThe")
Sim_OnsetThe_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfSim_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
Sim_OnsetThe_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType - Relatedness|participant) + (1+Relatedness|Verb), data = dfSim_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
Sim_OnsetThe_inverse <- glmer((OnsetTheRaw) ~ Relatedness*VerbType+ (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = subset(df, df$soa == 'pn0'), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) 

saveRDS(Sim_OnsetThe_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_log.rds")
saveRDS(Sim_OnsetThe_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_raw.rds")
saveRDS(Sim_OnsetThe_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_inverse.rds")

Sim_OnsetThe_log_coefs <- data.frame(coef(summary(Sim_OnsetThe_log)))
Sim_OnsetThe_log_coefs$p.z <- 2 * (1 - pnorm(abs(Sim_OnsetThe_log_coefs$t.value)))
Sim_OnsetThe_log_coefs
Sim_OnsetThe_raw_coefs <- data.frame(coef(summary(Sim_OnsetThe_raw)))
Sim_OnsetThe_raw_coefs$p.z <- 2 * (1 - pnorm(abs(Sim_OnsetThe_raw_coefs$t.value)))
Sim_OnsetThe_raw_coefs

#### comparing the interference effect between regions

dfNeg_unacc <- subset(dfNeg, dfNeg$VerbType == "Unacc")
Neg_OnsetThe_region_log <- lmer((ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region|participant) + (1+Relatedness*Region|Verb), data = dfNeg_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))
Neg_OnsetThe_region_raw <- lmer(exp(ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region|participant) + (1+Relatedness*Region|Verb), data = dfNeg_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))

df_regionRaw <- df %>% gather(Region, ProductionTime, OnsetTheRaw, dogRaw, aboveRaw, theRaw, appleRaw, isRaw ,AIRaw, ATRaw,postRaw,preNRaw)
df_regionRaw$Region <- factor(df_regionRaw$Region, ordered = T, levels = c("OnsetTheRaw", "dogRaw", "aboveRaw", "theRaw", "appleRaw", "isRaw", "AIRaw","ATRaw","postRaw","preNRaw"))
df_ROIRaw <- subset(df_regionRaw, df_regionRaw$Region == 'OnsetTheRaw'|df_regionRaw$Region == "postRaw")
df_ROIRaw$Region <- factor(df_ROIRaw$Region, ordered = T, levels = c("OnsetTheRaw", "postRaw")); contrasts(df_ROIRaw$Region) <- c(-0.5,0.5)
dfNeg_unaccRaw <- subset(df_ROIRaw, df_ROIRaw$soa == 'n150' & df_ROIRaw$VerbType == "Unacc")

Neg_OnsetThe_region_inverse <- glmer((ProductionTime) ~ Relatedness*Region + (1+Region+Relatedness|participant) + (1 + Relatedness + Region|Verb), data = dfNeg_unaccRaw, family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) ### including interaction terms in the random effect caused the convergence failure

saveRDS(Neg_OnsetThe_region_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_region_log.rds")
saveRDS(Neg_OnsetThe_region_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_region_raw.rds")
saveRDS(Neg_OnsetThe_region_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_region_inverse.rds")
Neg_OnsetThe_region_log  <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_region_log.rds")
Neg_OnsetThe_region_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_region_raw.rds")
Neg_OnsetThe_region_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Neg_OnsetThe_region_inverse.rds")


dfSim_unacc <- subset(dfSim, dfSim$VerbType == "Unacc")
Sim_OnsetThe_region_log <- lmer((ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region - Relatedness|participant) + (1+Relatedness*Region|Verb), data = dfSim_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))
Sim_OnsetThe_region_raw <- lmer(exp(ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region|participant) + (1+Relatedness*Region|Verb), data = dfSim_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))

dfSim_unaccRaw <- subset(df_ROIRaw, df_ROIRaw$soa == 'pn0' & df_ROIRaw$VerbType == "Unacc")
Sim_OnsetThe_region_inverse <- glmer((ProductionTime) ~ Relatedness*Region + (1+Region*Relatedness|participant) + (1 + Relatedness*Region|Verb), data = dfSim_unaccRaw, family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) ### including interaction terms in the random effect caused the convergence failure

saveRDS(Sim_OnsetThe_region_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_region_log.rds")
saveRDS(Sim_OnsetThe_region_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_region_raw.rds")
saveRDS(Sim_OnsetThe_region_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_region_inverse.rds")
Sim_OnsetThe_region_log  <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_region_log.rds")
Sim_OnsetThe_region_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_region_raw.rds")
Sim_OnsetThe_region_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment3/rmodels/Sim_OnsetThe_region_inverse.rds")



##### Second ROI

df_ROI <- subset(df_region, df_region$Region == 'pre'|df_region$Region == "AI")
df_ROI$Region <- factor(df_ROI$Region, ordered = T, levels = c("pre", "AI")); contrasts(df_ROI$Region) <- c(-0.5,0.5)

dfNeg <- subset(df_ROI, df_ROI$soa == "n150")
dfSim <- subset(df_ROI, df_ROI$soa == 'pn0')

dfNeg_AI <- subset(dfNeg, dfNeg$Region == "AI")
dfSim_AI <- subset(dfSim, dfSim$Region == "AI")

Neg_AI_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfNeg_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
Neg_AI_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfNeg_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
Neg_AI_inverse <- glmer((AIRaw) ~ Relatedness*VerbType + (1+ Relatedness*VerbType - Relatedness|participant) + (1+Relatedness|Verb), data = subset(df, df$soa == "n150"), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) #### including by-subject random slope of Relatedness caused convergence failure

Neg_AI_log_coefs <- data.frame(coef(summary(Neg_AI_log)))
Neg_AI_log_coefs$p.z <- 2 * (1 - pnorm(abs(Neg_AI_log_coefs$t.value)))
Neg_AI_log_coefs
Neg_AI_raw_coefs <- data.frame(coef(summary(Neg_AI_raw)))
Neg_AI_raw_coefs$p.z <- 2 * (1 - pnorm(abs(Neg_AI_raw_coefs$t.value)))
Neg_AI_raw_coefs

saveRDS(Neg_AI_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/rmodels/Neg_AI_log.rds")
saveRDS(Neg_AI_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/rmodels/Neg_AI_raw.rds")
saveRDS(Neg_AI_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/rmodels/Neg_AI_inverse.rds")


dfSim_AI <- subset(dfSim, dfSim$Region == "AI")

Sim_AI_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfSim_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
Sim_AI_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfSim_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
Sim_AI_inverse <- glmer((AIRaw) ~ Relatedness*VerbType + (1+ Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = subset(df, df$soa == "pn0"), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) #### including by-subject random slope of Relatedness caused convergence failure

Sim_AI_log_coefs <- data.frame(coef(summary(Sim_AI_log)))
Sim_AI_log_coefs$p.z <- 2 * (1 - pnorm(abs(Sim_AI_log_coefs$t.value)))
Sim_AI_log_coefs
Sim_AI_raw_coefs <- data.frame(coef(summary(Sim_AI_raw)))
Sim_AI_raw_coefs$p.z <- 2 * (1 - pnorm(abs(Sim_AI_raw_coefs$t.value)))
Sim_AI_raw_coefs

saveRDS(V_AI_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/rmodels/Sim_AI_log.rds")
saveRDS(V_AI_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/rmodels/Sim_AI_raw.rds")
saveRDS(V_AI_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/rmodels/Sim_AI_inverse.rds")

