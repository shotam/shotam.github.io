##### load data and packages #####

library('Rmisc'); library('lme4');library('multcomp');library('ggplot2');library('lmerTest');library('ordinal');library('tidyr');library('papaja')
df <- read.csv('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Codability/combined_data.csv')
df <- subset(df, df$experiment == "Exp3")
df$Picture<- interaction(df$Pic, df$Anoun)
df <- subset(df, df$participant != "211")
sum(is.na(df$Onset))/nrow(df)
df <- subset(df, df$Onset > 0)


error <- data.frame(table(df$Relatedness, df$VerbType, df$DistractorType))
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


#ggplot(subset(df, df$VerbType == "Unerg" & df$DistractorType == "Verb"), aes(AI, color=Relatedness)) +
#  geom_histogram(position="identity", binwidth=50, aes(y=..density.., fill=Relatedness),  alpha=0.5) +
#  geom_density()+
#  theme_classic()

##### coding factors #####

df$VerbType = factor(df$VerbType); contrasts(df$VerbType) <- c(-0.5,0.5)  ### unacc = -0.5 ; unerg = 0.5
df$Relatedness = factor(df$Relatedness); contrasts(df$Relatedness) <- c(-0.5,0.5) ### related = -0.5 ; unrelated = 0.5
df$DistractorType = factor(df$DistractorType); contrasts(df$DistractorType) <- c(-0.5,0.5) ### 0.5 = Verb, -0.5 = noun
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
df_region <- df %>% gather(Region, ProductionTime, Onset, The, dog, above, the, apple, is, OnsetThe ,AI, AT,post)
df_region$Region <- factor(df_region$Region, ordered = T, levels = c("Onset", "The" , "dog", "above", "the", "apple", "is", "OnsetThe" , "AI","AT","post"))

df_region_bySubj <- subset(df_region, df_region$Region != "OnsetThe" & df_region$Region != "AI" &df_region$Region != "AT" & df_region$Region != "post") %>% dplyr::group_by(participant,VerbType, Relatedness, DistractorType, Region) %>%
  dplyr::summarise(ProductionTime = mean(exp(ProductionTime), na.rm = TRUE))

descriptives <- df_region_bySubj %>% dplyr::group_by(DistractorType, VerbType, Relatedness, Region) %>%
  dplyr::summarise(
    Mean = as.integer(mean((ProductionTime), na.rm = TRUE)),
  )

descriptives <- spread(descriptives, Region, Mean)
apa_table(descriptives)

##### creating graph: region-by-region difference #####
dfV_graph <- subset(df, df$DistractorType == 'Verb')
graph1 <- as.data.frame(tapply(exp(dfV_graph$Onset),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph1$UnaccDiff <- graph1$Related.Unacc - graph1$Unrelated.Unacc
graph1$UnergDiff <- graph1$Related.Unerg - graph1$Unrelated.Unerg
graph1$region = '1'
graph2 <- as.data.frame(tapply(exp(dfV_graph$The),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph2$UnaccDiff <- graph2$Related.Unacc - graph2$Unrelated.Unacc
graph2$UnergDiff <- graph2$Related.Unerg - graph2$Unrelated.Unerg
graph2$region = '2'
graph3 <- as.data.frame(tapply(exp(dfV_graph$dog),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph3$UnaccDiff <- graph3$Related.Unacc - graph3$Unrelated.Unacc
graph3$UnergDiff <- graph3$Related.Unerg - graph3$Unrelated.Unerg
graph3$region = '3'
graph4 <- as.data.frame(tapply(exp(dfV_graph$above),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph4$UnaccDiff <- graph4$Related.Unacc - graph4$Unrelated.Unacc
graph4$UnergDiff <- graph4$Related.Unerg - graph4$Unrelated.Unerg
graph4$region = '4'
graph5 <- as.data.frame(tapply(exp(dfV_graph$the),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph5$UnaccDiff <- graph5$Related.Unacc - graph5$Unrelated.Unacc
graph5$UnergDiff <- graph5$Related.Unerg - graph5$Unrelated.Unerg
graph5$region = '5'
graph6 <- as.data.frame(tapply(exp(dfV_graph$apple),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph6$UnaccDiff <- graph6$Related.Unacc - graph6$Unrelated.Unacc
graph6$UnergDiff <- graph6$Related.Unerg - graph6$Unrelated.Unerg
graph6$region = '6'
graph7 <- as.data.frame(tapply(exp(dfV_graph$is),list(dfV_graph$participant,dfV_graph$Relatedness,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
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

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/VdifferenceXXX.pdf')
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


dfN_graph <- subset(df, df$DistractorType == 'Noun')
graph1 <- as.data.frame(tapply(exp(dfN_graph$Onset),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph1$UnaccDiff <- graph1$Related.Unacc - graph1$Unrelated.Unacc
graph1$UnergDiff <- graph1$Related.Unerg - graph1$Unrelated.Unerg
graph1$region = '1'
graph2 <- as.data.frame(tapply(exp(dfN_graph$The),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph2$UnaccDiff <- graph2$Related.Unacc - graph2$Unrelated.Unacc
graph2$UnergDiff <- graph2$Related.Unerg - graph2$Unrelated.Unerg
graph2$region = '2'
graph3 <- as.data.frame(tapply(exp(dfN_graph$dog),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph3$UnaccDiff <- graph3$Related.Unacc - graph3$Unrelated.Unacc
graph3$UnergDiff <- graph3$Related.Unerg - graph3$Unrelated.Unerg
graph3$region = '3'
graph4 <- as.data.frame(tapply(exp(dfN_graph$above),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph4$UnaccDiff <- graph4$Related.Unacc - graph4$Unrelated.Unacc
graph4$UnergDiff <- graph4$Related.Unerg - graph4$Unrelated.Unerg
graph4$region = '4'
graph5 <- as.data.frame(tapply(exp(dfN_graph$the),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph5$UnaccDiff <- graph5$Related.Unacc - graph5$Unrelated.Unacc
graph5$UnergDiff <- graph5$Related.Unerg - graph5$Unrelated.Unerg
graph5$region = '5'
graph6 <- as.data.frame(tapply(exp(dfN_graph$apple),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph6$UnaccDiff <- graph6$Related.Unacc - graph6$Unrelated.Unacc
graph6$UnergDiff <- graph6$Related.Unerg - graph6$Unrelated.Unerg
graph6$region = '6'
graph7 <- as.data.frame(tapply(exp(dfN_graph$is),list(dfN_graph$participant,dfN_graph$Relatedness,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
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

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Ndifference.pdf')
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
        axis.text.x  = element_text(angle=45, hjust = 1, size=14,face=c("bold","bold","italic","bold","bold","italic","italic"))) +
  coord_cartesian(ylim = c(-50, 100))
dev.off()

##### creating graphs: ROI analysis #####
library(Hmisc)
dfV_graph <- subset(dfV_graph, dfV_graph$OnsetThe > 0)
pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/VOnset.pdf')
apa_lineplot(
  data = dfV_graph,
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

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/VAI.pdf')
apa_lineplot(
  dfV_graph,
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

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/NOnset.pdf')
apa_lineplot(
  dfN_graph,
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

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/NAT.pdf')
apa_lineplot(
  dfN_graph,
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

##### VerbDistractor #####

df_ROI <- subset(df_region, df_region$Region == 'OnsetThe'|df_region$Region == "post")
df_ROI$Region <- factor(df_ROI$Region, ordered = T, levels = c("OnsetThe", "post")); contrasts(df_ROI$Region) <- c(-0.5,0.5)

dfV <- subset(df_ROI, df_ROI$DistractorType == 'Verb')

##### First ROI

dfV_OnsetThe <- subset(dfV, dfV$Region == "OnsetThe")
V_OnsetThe_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfV_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
V_OnsetThe_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType - Relatedness|participant) + (1+Relatedness|Verb), data = dfV_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
V_OnsetThe_inverse <- glmer((OnsetTheRaw) ~ Relatedness*VerbType + (1+Relatedness + VerbType|participant) + (1 + Relatedness|Verb), data = subset(df, df$DistractorType == "Verb"), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) #### including by-subject random slope of Relatedness caused convergence failure


saveRDS(V_OnsetThe_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_OnsetThe_log.rds")
saveRDS(V_OnsetThe_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_OnsetThe_raw.rds")
saveRDS(V_OnsetThe_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_OnsetThe_inverse.rds")
V_OnsetThe_log <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_OnsetThe_log.rds")
V_OnsetThe_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_OnsetThe_raw.rds")
V_OnsetThe_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_OnsetThe_inverse.rds")

V_OnsetThe_log_coefs <- data.frame(coef(summary(V_OnsetThe_log)))
V_OnsetThe_log_coefs$p.z <- 2 * (1 - pnorm(abs(V_OnsetThe_log_coefs$t.value)))
V_OnsetThe_log_coefs
V_OnsetThe_raw_coefs <- data.frame(coef(summary(V_OnsetThe_raw)))
V_OnsetThe_raw_coefs$p.z <- 2 * (1 - pnorm(abs(V_OnsetThe_raw_coefs$t.value)))
V_OnsetThe_raw_coefs


dfV_unacc <- subset(dfV, dfV$VerbType == "Unacc")
V_unacc_log <- lmer((ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region - Relatedness|participant) + (1 + Relatedness*Region|Verb), data = dfV_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))
V_unacc_raw <- lmer(exp(ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region|participant) + (1 + Relatedness*Region|Verb), data = dfV_unacc, control=lmerControl(optCtrl=list(maxfun=100000)))

V_unacc_log_coefs <- data.frame(coef(summary(V_unacc_log)))
V_unacc_log_coefs$p.z <- 2 * (1 - pnorm(abs(V_unacc_log_coefs$t.value)))
V_unacc_log_coefs
V_unacc_raw_coefs <- data.frame(coef(summary(V_unacc_raw)))
V_unacc_raw_coefs$p.z <- 2 * (1 - pnorm(abs(V_unacc_raw_coefs$t.value)))
V_unacc_raw_coefs

df_regionRaw <- df %>% gather(Region, ProductionTime, OnsetTheRaw, dogRaw, aboveRaw, theRaw, appleRaw, isRaw ,AIRaw, ATRaw,postRaw)
df_regionRaw$Region <- factor(df_regionRaw$Region, ordered = T, levels = c("OnsetTheRaw", "dogRaw", "aboveRaw", "theRaw", "appleRaw", "isRaw", "AIRaw","ATRaw","postRaw"))
df_ROIRaw <- subset(df_regionRaw, df_regionRaw$Region == 'OnsetTheRaw'|df_regionRaw$Region == "postRaw")
df_ROIRaw$Region <- factor(df_ROIRaw$Region, ordered = T, levels = c("OnsetTheRaw", "postRaw")); contrasts(df_ROIRaw$Region) <- c(-0.5,0.5)
dfVRaw <- subset(df_ROIRaw, df_ROIRaw$DistractorType == 'Verb')
dfV_unaccRaw <- subset(dfVRaw, dfVRaw$VerbType == "Unacc")
V_unacc_inverse <- glmer((ProductionTime) ~ Relatedness*Region + (1 + Relatedness*Region|participant) + (1 + Relatedness*Region|Verb), data = dfV_unaccRaw, family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa'))

saveRDS(V_unacc_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_unacc_log.rds")
saveRDS(V_unacc_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_unacc_raw.rds")
saveRDS(V_unacc_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_unacc_inverse.rds")

##### Second ROI

df_ROI <- subset(df_region, df_region$Region == 'pre'|df_region$Region == "AI")
df_ROI$Region <- factor(df_ROI$Region, ordered = T, levels = c("pre", "AI")); contrasts(df_ROI$Region) <- c(-0.5,0.5)
dfV <- subset(df_ROI, df_ROI$DistractorType == 'Verb')

dfV_AI <- subset(dfV, dfV$Region == "AI")
V_AI_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType - VerbType|participant) + (1+Relatedness|Verb), data = dfV_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
V_AI_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = dfV_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
V_AI_inverse <- glmer(AIRaw ~ Relatedness*VerbType + (1 + Relatedness|participant) + (1 + Relatedness|Verb), data = subset(df, df$DistractorType == "Verb"), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa')) #### including by-subject random slope of Relatedness caused convergence failure

saveRDS(V_AI_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_AI_log.rds")
saveRDS(V_AI_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_AI_raw.rds")
saveRDS(V_AI_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_AI_inverse.rds")
V_AI_log <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_AI_log.rds")
V_AI_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_AI_raw.rds")
V_AI_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/V_AI_inverse.rds")

V_AI_log_coefs <- data.frame(coef(summary(V_AI_log)))
V_AI_log_coefs$p.z <- 2 * (1 - pnorm(abs(V_AI_log_coefs$t.value)))
V_AI_log_coefs
V_AI_raw_coefs <- data.frame(coef(summary(V_AI_raw)))
V_AI_raw_coefs$p.z <- 2 * (1 - pnorm(abs(V_AI_raw_coefs$t.value)))
V_AI_raw_coefs


##### NounDistractor #####

df_region <- df %>% gather(Region, ProductionTime, OnsetThe, post)
df_region$Region <- factor(df_region$Region, levels = c("OnsetThe", "post"), ordered = TRUE)
contrasts(df_region$Region) <- c(-0.5, 0.5)
dfN <- subset(df_region, df_region$DistractorType == 'Noun') 
dfN_OnsetThe <- subset(dfN, dfN$Region == "OnsetThe")

N_OnsetThe_log <- lmer((ProductionTime) ~Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Verb), data = dfN_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
N_OnsetThe_raw <- lmer(exp(ProductionTime) ~Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Verb), data = dfN_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
N_OnsetThe_inverse <- glmer(OnsetTheRaw ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1+Relatedness|Verb), data = subset(df, df$DistractorType == "Noun"), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa'))

saveRDS(N_OnsetThe_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_OnsetThe_log.rds")
saveRDS(N_OnsetThe_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_OnsetThe_raw.rds")
saveRDS(N_OnsetThe_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_OnsetThe_inverse.rds")

N_OnsetThe_log <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_OnsetThe_log.rds")
N_OnsetThe_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_OnsetThe_raw.rds")
N_OnsetThe_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_OnsetThe_inverse.rds")


N_OnsetThe_log_coefs <- data.frame(coef(summary(N_OnsetThe_log)))
N_OnsetThe_log_coefs$p.z <- 2 * (1 - pnorm(abs(N_OnsetThe_log_coefs$t.value)))
N_OnsetThe_log_coefs
N_OnsetThe_raw_coefs <- data.frame(coef(summary(N_OnsetThe_raw)))
N_OnsetThe_raw_coefs$p.z <- 2 * (1 - pnorm(abs(N_OnsetThe_raw_coefs$t.value)))
N_OnsetThe_raw_coefs


df_region <- df %>% gather(Region, ProductionTime, AT, pre)
df_region$Region <- factor(df_region$Region, levels = c("pre", "AT"), ordered = TRUE)
contrasts(df_region$Region) <- c(-0.5, 0.5)
dfN <- subset(df_region, df_region$DistractorType == 'Noun') 

dfN_AT <- subset(dfN, dfN$Region == "AT")
N_AT_log <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Verb), data = dfN_AT, control=lmerControl(optCtrl=list(maxfun=100000)))
N_AT_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Verb), data = dfN_AT, control=lmerControl(optCtrl=list(maxfun=100000)))
N_AT_inverse <- glmer(ATRaw ~ Relatedness*VerbType + (1 + Relatedness*VerbType - Relatedness|participant) + (1+Relatedness|Verb), data = subset(df, df$DistractorType == "Noun"), family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa'))


saveRDS(N_AT_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_AT_log.rds")
saveRDS(N_AT_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_AT_raw.rds")
saveRDS(N_AT_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_AT_inverse.rds")

N_AT_log <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_AT_log.rds")
N_AT_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_AT_raw.rds")
N_AT_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_AT_inverse.rds")

N_AT_log_coefs <- data.frame(coef(summary(N_AT_log)))
N_AT_log_coefs$p.z <- 2 * (1 - pnorm(abs(N_AT_log_coefs$t.value)))
N_AT_log_coefs
N_AT_raw_coefs <- data.frame(coef(summary(N_AT_raw)))
N_AT_raw_coefs$p.z <- 2 * (1 - pnorm(abs(N_AT_raw_coefs$t.value)))
N_AT_raw_coefs


N_region_log <- lmer((ProductionTime) ~ Relatedness*VerbType*Region + (1 + Relatedness*VerbType*Region|participant) + (1 + Relatedness*Region|Verb), data = dfN, control=lmerControl(optCtrl=list(maxfun=100000)))
N_region_raw <- lmer(exp(ProductionTime) ~ Relatedness*VerbType*Region + (1 + Relatedness*VerbType*Region|participant) + (1 + Relatedness*Region|Verb), data = dfN, control=lmerControl(optCtrl=list(maxfun=100000)))

N_region_log_coefs <- data.frame(coef(summary(N_region_log)))
N_region_log_coefs$p.z <- 2 * (1 - pnorm(abs(N_region_log_coefs$t.value)))
N_region_log_coefs
N_region_raw_coefs <- data.frame(coef(summary(N_region_raw)))
N_region_raw_coefs$p.z <- 2 * (1 - pnorm(abs(N_region_raw_coefs$t.value)))
N_region_raw_coefs


df_regionRaw <- df %>% gather(Region, ProductionTime, OnsetTheRaw, dogRaw, aboveRaw, theRaw, appleRaw, isRaw ,AIRaw, ATRaw,postRaw,preNRaw)
df_regionRaw$Region <- factor(df_regionRaw$Region, ordered = T, levels = c("OnsetTheRaw", "dogRaw", "aboveRaw", "theRaw", "appleRaw", "isRaw", "AIRaw","ATRaw","postRaw","preNRaw"))
df_ROIRaw <- subset(df_regionRaw, df_regionRaw$Region == 'preNRaw'|df_regionRaw$Region == "ATRaw")
df_ROIRaw$Region <- factor(df_ROIRaw$Region, ordered = T, levels = c("preNRaw", "ATRaw")); contrasts(df_ROIRaw$Region) <- c(-0.5,0.5)
dfNRaw <- subset(df_ROIRaw, df_ROIRaw$DistractorType == 'Noun')

N_region_inverse <- glmer((ProductionTime) ~ Relatedness*VerbType*Region + (1 + Relatedness*Region|participant) + (1 + Relatedness*Region - Relatedness|Verb), data = dfNRaw, family = inverse.gaussian(link = "identity"), control=glmerControl(optCtrl=list(maxfun=100000), optimizer = 'bobyqa'))

saveRDS(N_region_log, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_region_log.rds")
saveRDS(N_region_raw, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_region_raw.rds")
saveRDS(N_region_inverse, file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_region_inverse.rds")

N_region_log <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_region_log.rds")
N_region_raw <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_region_raw.rds")
N_region_inverse <- readRDS(file = "/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/rmodels/N_region_inverse.rds")




