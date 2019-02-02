##### load data and packages #####
library('Rmisc'); library('lme4');library('multcomp');library('ggplot2');library('lmerTest');library('ordinal');library(tidyr)

df <- read.csv('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_data.csv', comment.char="#")
df$soa <- 'p300'
df <- subset(df, df$participant != '211')
df <- subset(df, df$participant != '284')
sum(is.na(df$Onset))/nrow(df)
df <- subset(df, df$Onset > 0)

##### exclude the extreme values #####
df$OnsetThe[df$OnsetThe > 5000 |df$OnsetThe < 300] <- NA ; sum(is.na(df$OnsetThe))/nrow(df)
df$Onset[df$Onset > 5000 |df$Onset < 300] <- NA; sum(is.na(df$Onset))/nrow(df)
df$The[df$The > 2000 |df$The < 30] <- NA; sum(is.na(df$The))/nrow(df)
df$dog[df$dog > 2000 |df$dog < 30 ] <- NA; sum(is.na(df$dog))/nrow(df)
df$above[df$above > 2000 |df$above < 30 ] <- NA; sum(is.na(df$above))/nrow(df)
df$the[df$the > 2000 |df$the < 30 ] <- NA;  sum(is.na(df$the))/nrow(df)
df$apple[df$apple > 2000 |df$apple < 30 ] <- NA;  sum(is.na(df$apple))/nrow(df)
df$is[df$is > 2000 |df$is < 30 ] <- NA;  sum(is.na(df$is))/nrow(df)
df$AI <- df$apple + df$is
df$AT <- df$above + df$the
#df <- subset(df, df$OnsetThe > 0 &  df$dog > 0 & df$above > 0 & df$the > 0 & df$apple >0 & df$is >0) 

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

##### coding factors #####

df$VerbType = factor(df$VerbType); contrasts(df$VerbType) <- c(-0.5,0.5)  ### unacc = -0.5 ; unerg = 0.5
df$Relatedness = factor(df$Relatedness); contrasts(df$Relatedness) <- c(-0.5,0.5) ### related = -0.5 ; unrelated = 0.5
df$DistractorType = factor(df$DistractorType); contrasts(df$DistractorType) <- c(-0.5,0.5) ### 0.5 = Verb, -0.5 = noun
df$participant = factor(df$participant)
df$Transcription = factor(df$Transcription)

##### descriptive stats #####
df_region <- df %>% gather(Region, ProductionTime, Onset, The, dog, above, the, apple, is, OnsetThe, AI, AT)
df_region$Region <- factor(df_region$Region, ordered = T, levels = c("Onset", "The", "dog", "above", "the", "apple", "is","OnsetThe", "AI","AT"))

df_region_bySubj <- subset(df_region, df_region$Region != "OnsetThe" & df_region$Region != "AI" &df_region$Region != "AT") %>% dplyr::group_by(participant,VerbType, Relatedness, DistractorType, Region) %>%
  dplyr::summarise(ProductionTime = mean(exp(ProductionTime), na.rm = TRUE))

descriptives <- df_region_bySubj %>% dplyr::group_by(DistractorType, VerbType, Relatedness, Region) %>%
  dplyr::summarise(
    Mean = as.integer(mean((ProductionTime), na.rm = TRUE)),
  )

descriptives <- spread(descriptives, Region, Mean)
apa_table(descriptives)

##### VerbDistractor #####

df_ROI <- subset(df_region, df_region$Region == 'OnsetThe'|df_region$Region == "AI")
df_ROI$Region <- factor(df_ROI$Region, ordered = T, levels = c("OnsetThe", "AI")); contrasts(df_ROI$Region) <- c(-0.5,0.5)

dfV <- subset(df_ROI, df_ROI$DistractorType == 'Verb')

V_model <- lmer((ProductionTime) ~ Relatedness*VerbType*Region + (1 + Relatedness*VerbType*Region|participant) + (1 + Relatedness*Region|Transcription), data = dfV, control=lmerControl(optCtrl=list(maxfun=100000)))
V_model_summary <- summary(V_model)
write.csv(V_model_summary$coefficients, '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_Vmodel_summary.csv')
save(V_model, file = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_Vmodel.rda')
load(file = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_Vmodel.rda')
apa_table(V_model_summary$coefficients)

dfV_OnsetThe <- subset(dfV, dfV$Region == "OnsetThe")
V_OnsetThe <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Transcription), data = dfV_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
V_OnsetThe_summary <- summary(V_OnsetThe)
apa_table(V_OnsetThe_summary$coefficients)
write.csv(V_OnsetThe_summary$coefficients, '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_V_OnsetThe_summary.csv')

dfV_AI <- subset(dfV, dfV$Region == "AI")
V_AI <- lmer(ProductionTime ~ Relatedness*VerbType+ (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Transcription), data = dfV_AI, control=lmerControl(optCtrl=list(maxfun=100000)))
V_AI_summary <- summary(V_AI)
apa_table(V_AI_summary$coefficients)
write.csv(V_AI_summary$coefficients, '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_V_AI_summary.csv')

##### Noun interference #####

df_region <- df %>% gather(Region, ProductionTime, OnsetThe, AT)
df_region$Region <- factor(df_region$Region, levels = c("OnsetThe", "AT"), ordered = TRUE)
contrasts(df_region$Region) <- c(-0.5, 0.5)
dfN <- subset(df_region, df_region$DistractorType == 'Noun') 
N_model <- lmer((ProductionTime) ~Relatedness*VerbType*Region + (1 + Relatedness*VerbType*Region|participant) + (1 + Relatedness*Region|Transcription), data = dfN, control=lmerControl(optCtrl=list(maxfun=100000)))
save(N_model, file = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_Nmodel.rda')
load(file = '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment1/Experiment2b_Vmodel.rda')
N_model_summary <- summary(N_model)
write.csv(N_model_summary$coefficients, '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_Nmodel_summary.csv')
apa_table(N_model_summary$coefficients)

dfN_OnsetThe <- subset(dfN, dfN$Region == "OnsetThe")
N_OnsetThe <- lmer((ProductionTime) ~Relatedness*VerbType + (1 + Relatedness+VerbType|participant) + (1 + Relatedness|Transcription), data = dfN_OnsetThe, control=lmerControl(optCtrl=list(maxfun=100000)))
N_OnsetThe_summary <- summary(N_OnsetThe)
apa_table(N_OnsetThe_summary$coefficients)
write.csv(N_OnsetThe_summary$coefficients, '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_N_OnsetThe_summary.csv')

dfN_AT <- subset(dfN, dfN$Region == "AT")
N_AT <- lmer((ProductionTime) ~ Relatedness*VerbType + (1 + Relatedness*VerbType|participant) + (1 + Relatedness|Transcription), data = dfN_AT, control=lmerControl(optCtrl=list(maxfun=100000)))
N_AT_summary <- summary(N_AT)
apa_table(N_AT_summary$coefficients)
write.csv(N_AT_summary$coefficients, '/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Experiment2b_N_AT_summary.csv')

##### creating graphs: ROI analysis #####

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/VOnset.pdf')
apa_lineplot(
  subset(dfV, (dfV$Region == 'OnsetThe')),
  id = "participant",
  dv = "ProductionTime",
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
  subset(dfV, (dfV$Region == 'AI')),
  id = "participant",
  dv = "ProductionTime",
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
  subset(dfN, (dfN$Region == 'OnsetThe')),
  id = "participant",
  dv = "ProductionTime",
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
  subset(dfN, (dfN$Region == 'AT')),
  id = "participant",
  dv = "ProductionTime",
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

##### creating graph: region-by-region interference #####
dfV_graph <- subset(df, df$DistractorType == 'Verb')
graph1 <- as.data.frame(tapply(exp(dfV_graph$Onset),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph1$UnaccDiff <- graph1$Vrel.Unacc - graph1$Vunrel.Unacc
graph1$UnergDiff <- graph1$Vrel.Unerg - graph1$Vunrel.Unerg
graph1$region = '1'
graph2 <- as.data.frame(tapply(exp(dfV_graph$The),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph2$UnaccDiff <- graph2$Vrel.Unacc - graph2$Vunrel.Unacc
graph2$UnergDiff <- graph2$Vrel.Unerg - graph2$Vunrel.Unerg
graph2$region = '2'
graph3 <- as.data.frame(tapply(exp(dfV_graph$dog),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph3$UnaccDiff <- graph3$Vrel.Unacc - graph3$Vunrel.Unacc
graph3$UnergDiff <- graph3$Vrel.Unerg - graph3$Vunrel.Unerg
graph3$region = '3'
graph4 <- as.data.frame(tapply(exp(dfV_graph$above),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph4$UnaccDiff <- graph4$Vrel.Unacc - graph4$Vunrel.Unacc
graph4$UnergDiff <- graph4$Vrel.Unerg - graph4$Vunrel.Unerg
graph4$region = '4'
graph5 <- as.data.frame(tapply(exp(dfV_graph$the),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph5$UnaccDiff <- graph5$Vrel.Unacc - graph5$Vunrel.Unacc
graph5$UnergDiff <- graph5$Vrel.Unerg - graph5$Vunrel.Unerg
graph5$region = '5'
graph6 <- as.data.frame(tapply(exp(dfV_graph$apple),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph6$UnaccDiff <- graph6$Vrel.Unacc - graph6$Vunrel.Unacc
graph6$UnergDiff <- graph6$Vrel.Unerg - graph6$Vunrel.Unerg
graph6$region = '6'
graph7 <- as.data.frame(tapply(exp(dfV_graph$is),list(dfV_graph$participant,dfV_graph$Condition,dfV_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph7$UnaccDiff <- graph7$Vrel.Unacc - graph7$Vunrel.Unacc
graph7$UnergDiff <- graph7$Vrel.Unerg - graph7$Vunrel.Unerg
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

pdf('/Users/shotamomma/Documents/GitHub/shotam.github.io/Experiments/Octopus/Experiment2b/Vdifference.pdf')
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
graph1 <- as.data.frame(tapply(exp(dfN_graph$Onset),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph1$UnaccDiff <- graph1$Nrel.Unacc - graph1$Nunrel.Unacc
graph1$UnergDiff <- graph1$Nrel.Unerg - graph1$Nunrel.Unerg
graph1$region = '1'
graph2 <- as.data.frame(tapply(exp(dfN_graph$The),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph2$UnaccDiff <- graph2$Nrel.Unacc - graph2$Nunrel.Unacc
graph2$UnergDiff <- graph2$Nrel.Unerg - graph2$Nunrel.Unerg
graph2$region = '2'
graph3 <- as.data.frame(tapply(exp(dfN_graph$dog),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph3$UnaccDiff <- graph3$Nrel.Unacc - graph3$Nunrel.Unacc
graph3$UnergDiff <- graph3$Nrel.Unerg - graph3$Nunrel.Unerg
graph3$region = '3'
graph4 <- as.data.frame(tapply(exp(dfN_graph$above),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph4$UnaccDiff <- graph4$Nrel.Unacc - graph4$Nunrel.Unacc
graph4$UnergDiff <- graph4$Nrel.Unerg - graph4$Nunrel.Unerg
graph4$region = '4'
graph5 <- as.data.frame(tapply(exp(dfN_graph$the),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph5$UnaccDiff <- graph5$Nrel.Unacc - graph5$Nunrel.Unacc
graph5$UnergDiff <- graph5$Nrel.Unerg - graph5$Nunrel.Unerg
graph5$region = '5'
graph6 <- as.data.frame(tapply(exp(dfN_graph$apple),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph6$UnaccDiff <- graph6$Nrel.Unacc - graph6$Nunrel.Unacc
graph6$UnergDiff <- graph6$Nrel.Unerg - graph6$Nunrel.Unerg
graph6$region = '6'
graph7 <- as.data.frame(tapply(exp(dfN_graph$is),list(dfN_graph$participant,dfN_graph$Condition,dfN_graph$VerbType), FUN = function(x) mean(x, na.rm = T)))
graph7$UnaccDiff <- graph7$Nrel.Unacc - graph7$Nunrel.Unacc
graph7$UnergDiff <- graph7$Nrel.Unerg - graph7$Nunrel.Unerg
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






