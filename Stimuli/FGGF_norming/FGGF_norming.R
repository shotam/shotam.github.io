#### change wd to your working directory
wd <- "~/Documents/Projects/FGGF_norming/"

df <- read.csv( paste(wd, "norming_data.csv", sep = ''))
df$match <- 0
df$match[as.character(df$picture) == as.character(df$verb_root)] <- 1

df$agreement <- ave(df$match, df$picture, FUN = function(x) sum(x, na.rm = T))
verbs <- levels(df$picture)

count = 0
mat = matrix()
mat2 = matrix()
for(i in verbs){
  count = count + 1
  sub <- subset(df, df$picture == i)
  sub$verb_root <- factor(sub$verb_root)
  mat[count] = c(sum(table(sub$verb_root)/length(sub$verb_root) * log2(1/(table(sub$verb_root)/length(sub$verb_root)))))
  mat2[count] = i
  codability = as.data.frame(cbind(mat2,mat))
}
colnames(codability) <- c('picture', 'Hstats')

agr <- aggregate(df$match, by = list(df$picture),FUN = function(x) sum(x, na.rm = T)/length(x))
colnames(agr) <- c('picture', 'agreement')

agr$Hstats <- codability$Hstats

write.csv(agr, file = paste(wd, "codability_norm.csv",sep = ''), row.names = FALSE)

