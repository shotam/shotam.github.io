recode.vars <- function(df, varList, codingFunction, scaleLevelDiffsTo1=F, maxNumLevels=10) {
  ###
  ###   df:                 REQUIRED argument: The data frame that contains the data to be submitted to the mixed effects model.
  ###
  ###   varList:            REQUIRED argument: A list of column names in the data frame that contain the variables you want to include in the mixed-effects model.
  ###
  ###   codingFunction:     REQUIRED argument: A function that will take an integer and return a matrix of variables to use for recoding.
  ###                                          Examples: contr.sum, contr.treatment, contr.helmert.
  ###
  ###   scaleLevelDiffsTo1: OPTIONAL argument: If you want the two relevant factor levels in each column to be separated by 1, set this argument to True.
  ###
  ###   maxNumLevels:       OPTIONAL argument: Don't worry about it unless varList includes a nominal variable with >10 levels; in that case, identify the
  ###                                          factor in varList that contains the most levels, and set this equal to (at least) that number of levels.
  ###
  
  # function to center a continuous variable
  ct <- function(x) {(scale(as.numeric(x), scale=F))}
  
  for (nextVar in varList) {
    
    # check to make sure variable is either a factor or is continuous (what else would it be?)
    if(!is.factor(df[[nextVar]]) & !is.numeric(df[[nextVar]])) { stop("Factor ", nextVar, " is not a factor and is not continuous.") }
    
    numLevels = length(names(summary(df[[nextVar]])))
    
    # recode if it's a factor variable with 3+ levels
    if(is.factor(df[[nextVar]])) {
      
      # don't execute if there are a ridiculous number of factor levels (unless explicitly passed a parameter allowing that many)
      #if(is.numeric(maxNumLevels)) { maxNumLevelsCheck = maxNumLevels } else { maxNumLevelsCheck = 10 }
      if(numLevels > maxNumLevels) {
        stop(paste("Factor ", nextVar, " has more levels (", numLevels, ") than currently permitted (", maxNumLevels, "). If this is intended, run recode.vars() with maxNumLevels >= ", numLevels, ".", sep=""))
      }
      
      # recode variable using supplied codingFunction
      recoded.matrix = sapply(df[[nextVar]], function(i) codingFunction( numLevels )[i,])
      
      # add each row to the df
      for (i in seq(numLevels-1)) {
        
        if(numLevels > 2) {
          colToAdd    <- recoded.matrix[i,]
          nextColName <- paste(nextVar, ".lev", i, sep="")
        } else {
          colToAdd    <- recoded.matrix
          nextColName <- paste(nextVar, ".lev",    sep="")
        }
        
        # check whether level differences should be scaled & scale if so
        if(scaleLevelDiffsTo1==T) {# & all(c(-1,1) %in% as.numeric(names(summary(as.factor(as.character(colToAdd))))))) {
          scalingFactor <- max(colToAdd) - min(colToAdd)
          colToAdd <- colToAdd/scalingFactor
        }
        
        df[[nextColName]] <- colToAdd
      }
      
      # center if it's a numeric/continuous variable
    } else {
      
      df[[paste(nextVar, ".ct", sep="")]] <- ct(df[[nextVar]])
      
    }
  }
  df
}