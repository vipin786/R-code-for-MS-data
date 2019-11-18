#Set working directory
setwd("G:/BSBE Project Nov2018/Proteomics_11Jan18")
################################
library(readxl)
library(caret)
library(glmnet)
library(dplyr)
library(ggplot2)
comb.data <- read_xlsx("Combined data analysis_2_09012019.xlsx", sheet = 4)
write.csv(comb.data, file = "comb_data.csv", row.names = FALSE)
df <- read.csv("comb_data.csv", stringsAsFactors = FALSE)
df.row.names <- df[,1]
df <- as.data.frame(t(df[,-1]))
colnames(df) <- df.row.names
dfcopy <- df[2:nrow(df),]
df <- df[5:nrow(df),]
df$Condition <- as.character(df$Condition)
df <- as.data.frame(df)
df.x <- df[,2:ncol(df)]
df.X <- as.numeric(as.matrix(df.x))
dim(df.X) <- dim(df.x)
row.names(df.X) <- row.names(df.x)
colnames(df.X) <- colnames(df.x)
df.new <- cbind.data.frame(df.X, df$Condition)
df.new$Condition <- df.new$`df$Condition`
df.new$`df$Condition` <- NULL
df <- df.new
summary(df)
df$Condition <- relevel(df$Condition, ref = "HC")

#number of cases for each category in the variable 'Condition'
with(df, table(Condition))
# Condition
# HC   CB  CCB  CSA   DF NSFM NSVM   SA   SD   SF  SVM 
# 18    6    4    4    4    8   12    6   12   12   20 
###############################
library(nnet)
library(Hmisc)
library(VIM)
library(FactoMineR)
library(missMDA)
library(naniar)

dim(na.omit(df))
gg_miss_var(df)

res<-summary(aggr(df, sortVar=TRUE))$combinations
varnum_missing <- aggr(df, sortVar=TRUE)
write.csv(varnum_missing$missings, "variable_num_missings.csv")

matrixplot(df, sort = 2)
#vis_miss(df, sort_miss = TRUE)

marginplot(df[,c("Condition","O60476")])
marginplot(df[,c("Condition","Q9UGM5")])

#percentage of missing values in data
pct_miss(df)
# [1] 17.57251

library(dplyr)
MissingValues <- df %>%
  group_by(df$Condition) %>%
  miss_var_summary()

write.csv(MissingValues, "missingvalues_info.csv")

#################################
#### Code to impute data ########
df.X.row.names <- row.names(df.X)
df.X_t <- as.data.frame(t(df.X))
colnames(df.X_t) <- paste("ALOG2", df.X.row.names, sep = ".")
summary(df.X_t)

log2.names <- grep("ALOG2", names(df.X_t), value = TRUE)
log2.names
conditions <- c("ALOG2.HC", "ALOG2.CCB", "ALOG2.CSA", "ALOG2.CB", "ALOG2.SA", "ALOG2.NSFM", "ALOG2.NSVM", "ALOG2.SF", "ALOG2.SVM", "ALOG2.DF", "ALOG2.SD")
min_count <- c(rep(2,length(conditions)))
cond.names <- lapply(conditions, # Group column names by conditions
                     function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
cond.names

cond.filter = sapply(1:length(cond.names), function(i) {
  df2 = df.X_t[cond.names[[i]]]   # Extract columns of interest
  df2 = as.matrix(df2)   # Cast as matrix for the following command
  sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
  sums >= min_count[i]   # Calculates whether min_count requirement is met
})

df.X_t$KEEP <- apply(cond.filter, 1, any)
df.X_t$KEEP

## Data imputation function
impute_data = function(df_imp, width = 0.3, downshift = 1.8) {
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  LOG2.names = grep("ALOG2", names(df_imp), value = TRUE)
  impute.names = gsub("ALOG2", "impute", LOG2.names)
  
  # Create new columns indicating whether the values are imputed
  df_imp[impute.names] = lapply(LOG2.names, function(x) !is.finite(df_imp[, x]))
  
  # Imputation
  set.seed(1)
  df_imp[LOG2.names] = lapply(LOG2.names,
                              function(x) {
                                temp = df_imp[[x]]
                                temp[!is.finite(temp)] = NA
                                
                                temp.sd = width * sd(temp[df_imp$KEEP], na.rm = TRUE)   # shrink sd width
                                temp.mean = mean(temp[df_imp$KEEP], na.rm = TRUE) - 
                                  downshift * sd(temp[df_imp$KEEP], na.rm = TRUE)   # shift mean of imputed values
                                
                                n.missing = sum(is.na(temp))
                                temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                                return(temp)
                              })
  return(df_imp)
}


## Apply imputation
df.FNI = impute_data(df.X_t)
write.csv(df.FNI, "Imputed_dataset.csv")


library(ggplot2)
ggplot(data = df.FNI, 
       aes(x = df.FNI$ALOG2.HC4)) + 
  geom_histogram() + 
  facet_wrap(df.FNI$impute.HC4)

glimpse(df.FNI)

#Imputed X (with rows as proteins and columns as condition)
df.ImputedX_t <- df.FNI[,log2.names]
df.ImputedX_t.rownames <- row.names(df.ImputedX_t)
#Imputed X
df.ImputedX <- as.data.frame(t(df.ImputedX_t))
colnames(df.ImputedX) <- df.ImputedX_t.rownames
#Imputed Data
df.Imputed <- df.ImputedX
df.Imputed$Condition <- df$Condition
summary(df.Imputed)
glimpse(df.Imputed)
write.csv(df.Imputed, "ImputedDatasetXY.csv")

##############################
# for (i in 0:10){
#   LogReg_HCD <- logRegCV(df_HCM, colname.response = "Condition", positive.response = "Dengue", k = 10, stratify = T, balance = "up", alpha = i/10)
#   LogReg_HCD_Results <- lapply(c('lambda.min','coef','pred'), function(i){ aggregateLogRegCV(LogReg_HCD, i) })
#   write.csv(rbind(LogReg_HCD_Results[[1]], LogReg_HCD_Results[[3]], LogReg_HCD_Results[[2]]), paste0("LogReg_HCD_Results", i, ".csv"))  
# }

# n <- 10
# k <- 10
# BootResamples <- createResample(data, times = n, list = TRUE)
# TrainTestSplits <- lapply(1:n, function(outerFold){
#   trainSet <- data[BootResamples[[outerFold]],]
#   testSet <- data[-unique(BootResamples[[outerFold]]),]
#   list(train = trainSet, test = testSet)
# })
# TrainTestSplits[[1]]$train
# TrainTestSplits[[1]]$test
#############################
alpha.grid = seq(0,1,by=0.1)
lambda.grid = 10^seq(2,-2,length=100)
srchGrid = expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)
kfoldcv = trainControl(method = "cv", number = 10, savePredictions = TRUE)
ENlogReg <- function(data = data, k = 10){
  OuterNFolds <- createFolds(data$Condition, k = k)
  en_model <- lapply(1:10, function(outerfold){
    cat('\n<< Outer fold CV: ', outerfold, ' >>')
    trainData <- data[-OuterNFolds[[outerfold]],]
    testData <- data[OuterNFolds[[outerfold]],]
    cat(" train : ", nrow(trainData), "test : ", nrow(testData))
    
    ### inner cross validation and prediction
    en.train <- train(Condition~., data = trainData, method = "glmnet", tuneGrid = srchGrid, metric = "Kappa", trControl = kfoldcv, standardize = TRUE)
    en.besttune <- en.train$bestTune
    en.predict <- predict(en.train, newdata = testData, s = en.train$bestTune$lambda)
    en.varImp <- varImp(en.train, s = en.train$bestTune$lambda)
    en.coef <- coef(en.train$finalModel, s = en.train$bestTune$lambda)
    en.confusionMatrix <- confusionMatrix(en.predict, as.factor(testData$Condition))
    en_model_results <- list(train = trainData, test = testData, trainModel = en.train, bestTune = en.besttune, predict = en.predict, varImp = en.varImp, confusionMatrix = en.confusionMatrix, coef = en.coef)
  })
  return(en_model)
}

kfoldcv2 = trainControl(method = "cv", number = 5, savePredictions = TRUE)
ENlogReg2 <- function(data = data, k = 5){
  OuterNFolds <- createFolds(data$Condition, k = k)
  en_model <- lapply(1:5, function(outerfold){
    cat('\n<< Outer fold CV: ', outerfold, ' >>')
    trainData <- data[-OuterNFolds[[outerfold]],]
    testData <- data[OuterNFolds[[outerfold]],]
    cat(" train : ", nrow(trainData), "test : ", nrow(testData))
    
    ### inner cross validation and prediction
    en.train <- train(Condition~., data = trainData, method = "glmnet", tuneGrid = srchGrid, metric = "Kappa", trControl = kfoldcv2, standardize = TRUE)
    en.besttune <- en.train$bestTune
    en.predict <- predict(en.train, newdata = testData, s = en.train$bestTune$lambda)
    en.varImp <- varImp(en.train, s = en.train$bestTune$lambda)
    en.coef <- coef(en.train$finalModel, s = en.train$bestTune$lambda)
    en.confusionMatrix <- confusionMatrix(en.predict, as.factor(testData$Condition))
    en_model_results <- list(train = trainData, test = testData, trainModel = en.train, bestTune = en.besttune, predict = en.predict, varImp = en.varImp, confusionMatrix = en.confusionMatrix, coef = en.coef)
  })
  return(en_model)
}

kfoldcv3 = trainControl(method = "cv", number = 4, savePredictions = TRUE)
ENlogReg3 <- function(data = data, k = 4){
  OuterNFolds <- createFolds(data$Condition, k = k)
  en_model <- lapply(1:4, function(outerfold){
    cat('\n<< Outer fold CV: ', outerfold, ' >>')
    trainData <- data[-OuterNFolds[[outerfold]],]
    testData <- data[OuterNFolds[[outerfold]],]
    cat(" train : ", nrow(trainData), "test : ", nrow(testData))
    
    ### inner cross validation and prediction
    en.train <- train(Condition~., data = trainData, method = "glmnet", tuneGrid = srchGrid, metric = "Kappa", trControl = kfoldcv3, standardize = TRUE)
    en.besttune <- en.train$bestTune
    en.predict <- predict(en.train, newdata = testData, s = en.train$bestTune$lambda)
    en.varImp <- varImp(en.train, s = en.train$bestTune$lambda)
    en.coef <- coef(en.train$finalModel, s = en.train$bestTune$lambda)
    en.confusionMatrix <- confusionMatrix(en.predict, as.factor(testData$Condition))
    en_model_results <- list(train = trainData, test = testData, trainModel = en.train, bestTune = en.besttune, predict = en.predict, varImp = en.varImp, confusionMatrix = en.confusionMatrix, coef = en.coef)
  })
  return(en_model)
}

aggregateENModel <- function(ENlogReg.object, var){
  ## coef
  if(var == 'coef'){
    output <- lapply(ENlogReg.object, function(fold){ fold$coef }) %>% do.call(cbind,.) %>% as.matrix() %>% as.data.frame()
    colnames(output) <- 1:ncol(output)
    class(output) <- c(class(output),'coef')
  }
  
  ## variable importance
  if(var == 'varImp'){
    output <- lapply(ENlogReg.object, function(fold){ fold$varImp$importance }) %>% do.call(cbind,.) %>% as.matrix() %>% as.data.frame()
    colnames(output) <- 1:ncol(output)
    class(output) <- c(class(output),'varImp')
  }
  
  ## confusionMatrixbyClass
  else if(var == "cMbyClass"){
    output <- lapply(ENlogReg.object, function(fold){ fold$confusionMatrix$byClass }) %>% do.call(cbind,.) %>% as.matrix() %>% as.data.frame()
    colnames(output) <- 1:ncol(output)
    class(output) <- c(class(output),'cMbyClass')
  }
  
  ## kappa
  else if(var == "kappa"){
    output <- lapply(ENlogReg.object, function(fold){ fold$confusionMatrix$overall }) %>% do.call(cbind,.) %>% as.matrix() %>% as.data.frame()
    colnames(output) <- 1:ncol(output)
    class(output) <- c(class(output),'kappa')
  }
  
  ## NumTrain
  else if(var == 'NumTrain'){
    output <- sapply(ENlogReg.object, function(fold){ nrow(fold$train) }) %>% t() %>% as.matrix() %>% as.data.frame()
    rownames(output) <- "NumTrain"
    colnames(output) <- 1:ncol(output)
  }
  
  ## NumTest
  else if(var == 'NumTest'){
    output <- sapply(ENlogReg.object, function(fold){ nrow(fold$test) }) %>% t() %>% as.matrix() %>% as.data.frame()
    rownames(output) <- "NumTest"
    colnames(output) <- 1:ncol(output)
  }
  
  ## alpha
  else if(var == 'alpha'){
    output <- sapply(ENlogReg.object, function(fold){ fold$bestTune$alpha }) %>% t() %>% as.matrix() %>% as.data.frame()
    rownames(output) <- "alpha"
    colnames(output) <- 1:ncol(output)
  }
  
  ## lambda
  else if(var == 'lambda'){
    output <- sapply(ENlogReg.object, function(fold){ fold$bestTune$lambda }) %>% t() %>% as.matrix() %>% as.data.frame()
    rownames(output) <- "lambda"
    colnames(output) <- 1:ncol(output)
  }
  # as.matrix(with(ENlogReg_HCD[[1]]$test, table(Condition)))
  
  ## TrainC
  else if(var == 'TrainC'){
    output <- sapply(ENlogReg.object, function(fold){ with(fold$train,table(Condition)) }) %>% as.matrix() %>% as.data.frame()
    colnames(output) <- 1:ncol(output)
  }
  
  ## TestC
  else if(var == 'TestC'){
    output <- sapply(ENlogReg.object, function(fold){ with(fold$test,table(Condition)) })%>% as.matrix() %>% as.data.frame()
    colnames(output) <- 1:ncol(output)
  }
  
  return(output)
}

######################################
significant.proteins <- read_xlsx("SignificantProteins.xlsx", sheet = 2)
write.csv(significant.proteins, file = "significantProteins.csv", row.names = FALSE)
significantProteins <- read.csv("significantProteins.csv", stringsAsFactors = FALSE)

#########################################
### HC vs Dengue vs Malaria dataset ####
df <- df.Imputed
# Condition
# HC   CB  CCB  CSA   DF NSFM NSVM   SA   SD   SF  SVM 
# 18    6    4    4    4    8   12    6   12   12   20
df$Condition <- as.character(df$Condition)
df <- df[df$Condition != "CCB",]
df <- df[df$Condition != "CSA",]
df$Condition <- gsub("CB", "Malaria", df$Condition)
df$Condition <- gsub("NSFM", "Malaria", df$Condition)
df$Condition <- gsub("NSVM", "Malaria", df$Condition)
df$Condition <- gsub("SA", "Malaria", df$Condition)
df$Condition <- gsub("SF", "Malaria", df$Condition)
df$Condition <- gsub("SVM", "Malaria", df$Condition)
df$Condition <- gsub("SD", "Dengue", df$Condition)
df$Condition <- gsub("DF", "Dengue", df$Condition)
df$Condition <- as.factor(df$Condition)
df$Condition <- relevel(df$Condition, ref = "HC")
df$Condition
with(df, table(Condition))
# Condition
# HC  Dengue Malaria 
# 18      16      64

#####################################################
### HC vs Dengue Elastic Net Logistic Regression ####
df_HCD <- df[df$Condition != "Malaria",]
df_HCD$Condition <- as.character(df_HCD$Condition)
df_HCD$Condition <- as.factor(df_HCD$Condition)
df_HCD$Condition <- relevel(df_HCD$Condition, ref = "HC")
with(df_HCD, table(Condition))

significantProteins.HCD <- significantProteins[significantProteins$DF == "DF",]$Protein
significantProteins.HCD
df_HCD <- df_HCD[c(significantProteins.HCD, "Condition")]

ENlogReg_HCD <- ENlogReg(data = df_HCD, k = 10)
ENlogReg_HCD_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCD, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCD_Results, "ENlogReg_HCD_Results.csv")

library(pROC)
ENlogReg_CN <- colnames(df_HCD)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_HCD, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$B1AHL2, col = 'red')
lines(roc$P49747, col = 'blue')
lines(roc$P01031, col = 'green')
legend(x = "bottomright", legend = paste(row.names(auc)[1:3], " (AUC = ", round(auc$AUC[1:3],4), ")", sep = ""),
       col = c("red","blue","green"), lty = "solid", cex=0.7)

plot(roc$B1AHL2, col = 'red')
lines(roc$P49747, col = 'blue')
lines(roc$P01031, col = 'green')
lines(roc$Q9UGM5, col = 'yellow')
lines(roc$P10643, col = 'purple')
legend(x = "bottomright", legend = paste(row.names(auc)[1:5], " (AUC = ", round(auc$AUC[1:5],4), ")", sep = ""),
       col = c("red","blue","green","yellow","purple"), lty = "solid", cex=0.7)

# roc <- roc(Condition ~ P10643 + B1AHL2 + P01031, data = df_HCD, family = 'binomial')
# plot(roc$P10643, col = 'red')
# lines(roc$B1AHL2, col = 'blue')
# lines(roc$P01031, col = 'green')
# legend(x = "bottomright", legend=c(paste("B1AHL2", " (auc = ", round(roc$B1AHL2$auc, 4), ")", sep = ""), paste("P01031", " (auc = ", round(roc$P01031$auc, 4), ")", sep = ""), paste("P10643", " (auc = ", round(roc$P10643$auc,4), ")", sep = "")),
#        col=c("blue", "green", "red"), lty=1, cex=0.8)
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

#########################################
### HC vs Malaria Logistic Regression ####
df_HCM <- df[df$Condition != "Dengue",]
df_HCM$Condition <- as.character(df_HCM$Condition)
df_HCM$Condition <- as.factor(df_HCM$Condition)
df_HCM$Condition <- relevel(df_HCM$Condition, ref = "HC")
with(df_HCM, table(Condition))

significantProteins.HCM <- significantProteins[significantProteins$VM == "VM" | significantProteins$FM == "FM",]$Protein
significantProteins.HCM
df_HCM <- df_HCM[c(significantProteins.HCM, "Condition")]

ENlogReg_HCM <- ENlogReg(data = df_HCM, k = 10)
ENlogReg_HCM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCM_Results, "ENlogReg_HCM_Results.csv")

ENlogReg_CN <- colnames(df_HCM)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_HCM, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$V9GYM3, col = 'red')
lines(roc$P02647, col = 'blue')
lines(roc$B4E1Z4, col = 'green')
legend(x = "bottomright", legend = paste(row.names(auc)[1:3], " (AUC = ", round(auc$AUC[1:3],4), ")", sep = ""),
       col = c("red","blue","green"), lty = 1, cex=0.7)

plot(roc$V9GYM3, col = 'red')
lines(roc$P02647, col = 'blue')
lines(roc$B4E1Z4, col = 'green')
lines(roc$P00738, col = 'yellow')
lines(roc$P02741, col = 'purple')
legend(x = "bottomright", legend = paste(row.names(auc)[1:5], " (AUC = ", round(auc$AUC[1:5],4), ")", sep = ""),
       col = c("red","blue","green","yellow","purple"), lty = 1, cex=0.7)

#########################################
### Dengue vs Malaria Logistic Regression ####
df_DM <- df[df$Condition != "HC",]
df_DM$Condition <- as.character(df_DM$Condition)
df_DM$Condition <- as.factor(df_DM$Condition)
df_DM$Condition <- relevel(df_DM$Condition, ref = "Dengue")
with(df_DM, table(Condition))

ENlogReg_DM <- ENlogReg(data = df_DM, k = 10)
ENlogReg_DM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_DM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_DM_Results, "ENlogReg_DM_Results.csv")

library(pROC)
ENlogReg_CN <- colnames(df_DM)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_DM, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$P02750, col = 'red')
lines(roc$P04040, col = 'yellow')
lines(roc$P06681, col = 'green')
lines(roc$P00450, col = 'blue')
legend(x = "bottomright", legend=c(paste("P00450", " (AUC = ", round(roc$P00450$auc,4), ")", sep = ""),paste("P02750", " (AUC = ", round(roc$P02750$auc, 4), ")", sep = ""), paste("P04040", " (AUC = ", round(roc$P04040$auc, 4), ")", sep = ""), paste("P06681", " (AUC = ", round(roc$P06681$auc,4), ")", sep = "")),
       col=c("blue", "red", "yellow", "green"), lty=1, cex=0.75, lwd = 2)

# plot(roc$B1AHL2, col = 'red')
# lines(roc$P49747, col = 'blue')
# lines(roc$P00740, col = 'green')
# lines(roc$C9JF17, col = 'yellow')
# lines(roc$P0DJI9, col = 'purple')
# legend(x = "bottomright", legend = paste(row.names(auc)[1:5], " (AUC = ", round(auc$AUC[1:5],4), ")", sep = ""),
#        col = c("red","blue","green","yellow","purple"), lty = 1, cex=0.7)
#########################################
#########################################
### Malaria - Falciparum vs Vivax #######
df <- df.Imputed
# Condition
# HC   CB  CCB  CSA   DF NSFM NSVM   SA   SD   SF  SVM 
# 18    6    4    4    4    8   12    6   12   12   20
df$Condition <- as.character(df$Condition)
df <- df[df$Condition != "CCB",]
df <- df[df$Condition != "CSA",]
df <- df[df$Condition != "SD",]
df <- df[df$Condition != "DF",]
df$Condition <- gsub("CB", "FM", df$Condition)
df$Condition <- gsub("NSFM", "FM", df$Condition)
df$Condition <- gsub("NSVM", "VM", df$Condition)
df$Condition <- gsub("SA", "FM", df$Condition)
df$Condition <- gsub("SF", "FM", df$Condition)
df$Condition <- gsub("SVM", "VM", df$Condition)
df$Condition <- as.factor(df$Condition)
df$Condition <- relevel(df$Condition, ref = "HC")
df$Condition
with(df, table(Condition))
# Condition
# HC FM VM 
# 18 32 32 

#####################################################
### HC vs FM Elastic Net Logistic Regression ####
df_HCFM <- df[df$Condition != "VM",]
df_HCFM$Condition <- as.character(df_HCFM$Condition)
df_HCFM$Condition <- as.factor(df_HCFM$Condition)
df_HCFM$Condition <- relevel(df_HCFM$Condition, ref = "HC")
with(df_HCFM, table(Condition))

significantProteins.HCFM <- significantProteins[significantProteins$FM == "FM",]$Protein
significantProteins.HCFM
df_HCFM <- df_HCFM[c(significantProteins.HCFM, "Condition")]

ENlogReg_HCFM <- ENlogReg(data = df_HCFM, k = 10)
ENlogReg_HCFM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCFM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCFM_Results, "ENlogReg_HCFM_Results.csv")

#####################################################
### HC vs VM Elastic Net Logistic Regression ####
df_HCVM <- df[df$Condition != "FM",]
df_HCVM$Condition <- as.character(df_HCVM$Condition)
df_HCVM$Condition <- as.factor(df_HCVM$Condition)
df_HCVM$Condition <- relevel(df_HCVM$Condition, ref = "HC")
with(df_HCVM, table(Condition))

significantProteins.HCVM <- significantProteins[significantProteins$VM == "VM",]$Protein
significantProteins.HCVM
df_HCVM <- df_HCVM[c(significantProteins.HCVM, "Condition")]

ENlogReg_HCVM <- ENlogReg(data = df_HCVM, k = 10)
ENlogReg_HCVM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCVM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCVM_Results, "ENlogReg_HCVM_Results.csv")

#####################################################
### FM vs VM Elastic Net Logistic Regression ####
df_VMFM <- df[df$Condition != "HC",]
df_VMFM$Condition <- as.character(df_VMFM$Condition)
df_VMFM$Condition <- as.factor(df_VMFM$Condition)
df_VMFM$Condition <- relevel(df_VMFM$Condition, ref = "VM")
with(df_VMFM, table(Condition))

significantProteins.VMFM <- significantProteins[significantProteins$VM == "VM" | significantProteins$FM == "FM",]$Protein
significantProteins.VMFM
df_VMFM <- df_VMFM[c(significantProteins.VMFM, "Condition")]

ENlogReg_VMFM <- ENlogReg(data = df_VMFM, k = 10)
ENlogReg_VMFM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_VMFM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_VMFM_Results, "ENlogReg_VMFM_Results.csv")

ENlogReg_CN <- colnames(df_VMFM)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_VMFM, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$P04040, col = 'red')
lines(roc$P01011, col = 'blue')
lines(roc$P25311, col = 'purple')
lines(roc$P04196, col = 'green')
lines(roc$P00450, col = 'yellow')
legend(x = "bottomright", legend=c(paste("P04040", " (AUC = ", round(roc$P04040$auc,4), ")", sep = ""),
                                   paste("P01011", " (AUC = ", round(roc$P01011$auc, 4), ")", sep = ""), 
                                   paste("P25311", " (AUC = ", round(roc$P25311$auc, 4), ")", sep = ""), 
                                   paste("P04196", " (AUC = ", round(roc$P04196$auc,4), ")", sep = ""),
                                   paste("P00450", " (AUC = ", round(roc$P00450$auc,4), ")", sep = "")),
       col=c("red", "blue", "purple", "green", "yellow"), lty=1, cex=0.75, lwd = 2)

###################################################################
###################################################################
### Falciparum Malaria Severity - HC vs NSFM vs CB vs SA vs SF ####
df <- df.Imputed
# Condition
# HC   CB  CCB  CSA   DF NSFM NSVM   SA   SD   SF  SVM 
# 18    6    4    4    4    8   12    6   12   12   20
df$Condition <- as.character(df$Condition)
df <- df[df$Condition == "HC" | df$Condition == "NSFM" | df$Condition == "CB" | df$Condition == "SA" | df$Condition == "SF",]
df$Condition <- as.factor(df$Condition)
df$Condition <- relevel(df$Condition, ref = "HC")
df$Condition
with(df, table(Condition))
# Condition
# HC  NSFM CB  SA   SF 
# 18  6    6   6   12 

#######################################################
### HC vs NSFM Elastic Net Logistic Regression ########
df_HCNSFM <- df[df$Condition == "HC" | df$Condition == "NSFM",]
df_HCNSFM$Condition <- as.character(df_HCNSFM$Condition)
df_HCNSFM$Condition <- as.factor(df_HCNSFM$Condition)
df_HCNSFM$Condition <- relevel(df_HCNSFM$Condition, ref = "HC")
with(df_HCNSFM, table(Condition))

significantProteins.HCNSFM <- significantProteins[significantProteins$NSFM == "NSFM",]$Protein
significantProteins.HCNSFM
df_HCNSFM <- df_HCNSFM[c(significantProteins.HCNSFM, "Condition")]

ENlogReg_HCNSFM <- ENlogReg2(data = df_HCNSFM, k = 5)
ENlogReg_HCNSFM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCNSFM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCNSFM_Results, "ENlogReg_HCNSFM_Results.csv")

#####################################################
### HC vs CB Elastic Net Logistic Regression ########
df_HCCB <- df[df$Condition == "HC" | df$Condition == "CB",]
df_HCCB$Condition <- as.character(df_HCCB$Condition)
df_HCCB$Condition <- as.factor(df_HCCB$Condition)
df_HCCB$Condition <- relevel(df_HCCB$Condition, ref = "HC")
with(df_HCCB, table(Condition))

significantProteins.HCCB <- significantProteins[significantProteins$CB == "CB",]$Protein
significantProteins.HCCB
df_HCCB <- df_HCCB[c(significantProteins.HCCB, "Condition")]

ENlogReg_HCCB <- ENlogReg2(data = df_HCCB, k = 5)
ENlogReg_HCCB_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCCB, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCCB_Results, "ENlogReg_HCCB_Results.csv")

#####################################################
### HC vs SA Elastic Net Logistic Regression ########
df_HCSA <- df[df$Condition == "HC" | df$Condition == "SA",]
df_HCSA$Condition <- as.character(df_HCSA$Condition)
df_HCSA$Condition <- as.factor(df_HCSA$Condition)
df_HCSA$Condition <- relevel(df_HCSA$Condition, ref = "HC")
with(df_HCSA, table(Condition))

significantProteins.HCSA <- significantProteins[significantProteins$SA == "SA",]$Protein
significantProteins.HCSA
df_HCSA <- df_HCSA[c(significantProteins.HCSA, "Condition")]

ENlogReg_HCSA <- ENlogReg2(data = df_HCSA, k = 5)
ENlogReg_HCSA_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCSA, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCSA_Results, "ENlogReg_HCSA_Results.csv")

#####################################################
### HC vs SF Elastic Net Logistic Regression ########
df_HCSF <- df[df$Condition == "HC" | df$Condition == "SF",]
df_HCSF$Condition <- as.character(df_HCSF$Condition)
df_HCSF$Condition <- as.factor(df_HCSF$Condition)
df_HCSF$Condition <- relevel(df_HCSF$Condition, ref = "HC")
with(df_HCSF, table(Condition))

significantProteins.HCSF <- significantProteins[significantProteins$SF == "SF",]$Protein
significantProteins.HCSF
df_HCSF <- df_HCSF[c(significantProteins.HCSF, "Condition")]

ENlogReg_HCSF <- ENlogReg(data = df_HCSF, k = 10)
ENlogReg_HCSF_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCSF, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCSF_Results, "ENlogReg_HCSF_Results.csv")

#####################################################
### NSFM vs CB Elastic Net Logistic Regression ########
df_NSFMCB <- df[df$Condition == "NSFM" | df$Condition == "CB",]
df_NSFMCB$Condition <- as.character(df_NSFMCB$Condition)
df_NSFMCB$Condition <- as.factor(df_NSFMCB$Condition)
df_NSFMCB$Condition <- relevel(df_NSFMCB$Condition, ref = "NSFM")
with(df_NSFMCB, table(Condition))

significantProteins.NSFMCB <- significantProteins[significantProteins$NSFM == "NSFM" | significantProteins$CB == "CB",]$Protein
significantProteins.NSFMCB
df_NSFMCB <- df_NSFMCB[c(significantProteins.NSFMCB, "Condition")]

ENlogReg_NSFMCB <- ENlogReg2(data = df_NSFMCB, k = 5)
ENlogReg_NSFMCB_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_NSFMCB, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_NSFMCB_Results, "ENlogReg_NSFMCB_Results.csv")

ENlogReg_CN <- colnames(df_NSFMCB)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_NSFMCB, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$P02647, col = 'red')
legend(x = "bottomright", legend=c(paste("P02647", " (AUC = ", round(roc$P02647$auc,4), ")", sep = "")),
       col=c("red"), lty=1, cex=0.75, lwd = 2)


#####################################################
### NSFM vs SA Elastic Net Logistic Regression ########
df_NSFMSA <- df[df$Condition == "NSFM" | df$Condition == "SA",]
df_NSFMSA$Condition <- as.character(df_NSFMSA$Condition)
df_NSFMSA$Condition <- as.factor(df_NSFMSA$Condition)
df_NSFMSA$Condition <- relevel(df_NSFMSA$Condition, ref = "NSFM")
with(df_NSFMSA, table(Condition))

significantProteins.NSFMSA <- significantProteins[significantProteins$NSFM == "NSFM" | significantProteins$SA == "SA",]$Protein
significantProteins.NSFMSA
df_NSFMSA <- df_NSFMSA[c(significantProteins.NSFMSA, "Condition")]

ENlogReg_NSFMSA <- ENlogReg2(data = df_NSFMSA, k = 5)
ENlogReg_NSFMSA_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_NSFMSA, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_NSFMSA_Results, "ENlogReg_NSFMSA_Results.csv")

ENlogReg_CN <- colnames(df_NSFMSA)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_NSFMSA, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$P02671, col = 'red')
lines(roc$P04040, col = 'blue')
legend(x = "bottomright", legend=c(paste("P02671", " (AUC = ", round(roc$P02671$auc,4), ")", sep = ""),
                                   paste("P04040", " (AUC = ", round(roc$P04040$auc, 4), ")", sep = "")),
       col=c("red", "blue"), lty=1, cex=0.75, lwd = 2)

#####################################################
### NSFM vs SF Elastic Net Logistic Regression ########
df_NSFMSF <- df[df$Condition == "NSFM" | df$Condition == "SF",]
df_NSFMSF$Condition <- as.character(df_NSFMSF$Condition)
df_NSFMSF$Condition <- as.factor(df_NSFMSF$Condition)
df_NSFMSF$Condition <- relevel(df_NSFMSF$Condition, ref = "NSFM")
with(df_NSFMSF, table(Condition))

significantProteins.NSFMSF <- significantProteins[significantProteins$NSFM == "NSFM" | significantProteins$SF == "SF",]$Protein
significantProteins.NSFMSF
df_NSFMSF <- df_NSFMSF[c(significantProteins.NSFMSF, "Condition")]

ENlogReg_NSFMSF <- ENlogReg2(data = df_NSFMSF, k = 5)
ENlogReg_NSFMSF_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_NSFMSF, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_NSFMSF_Results, "ENlogReg_NSFMSF_Results.csv")

ENlogReg_CN <- colnames(df_NSFMSF)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_NSFMSF, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$P02671, col = 'red')
lines(roc$P25311, col = 'blue')
legend(x = "bottomright", legend=c(paste("P02671", " (AUC = ", round(roc$P02671$auc,4), ")", sep = ""),
                                   paste("P25311", " (AUC = ", round(roc$P25311$auc, 4), ")", sep = "")),
       col=c("red", "blue"), lty=1, cex=0.75, lwd = 2)


#####################################################
### CB vs SA Elastic Net Logistic Regression ########
df_CBSA <- df[df$Condition == "CB" | df$Condition == "SA",]
df_CBSA$Condition <- as.character(df_CBSA$Condition)
df_CBSA$Condition <- as.factor(df_CBSA$Condition)
df_CBSA$Condition <- relevel(df_CBSA$Condition, ref = "CB")
with(df_CBSA, table(Condition))

significantProteins.CBSA <- significantProteins[significantProteins$CB == "CB" | significantProteins$SA == "SA",]$Protein
significantProteins.CBSA
df_CBSA <- df_CBSA[c(significantProteins.CBSA, "Condition")]

ENlogReg_CBSA <- ENlogReg2(data = df_CBSA, k = 5)
ENlogReg_CBSA_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_CBSA, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_CBSA_Results, "ENlogReg_CBSA_Results.csv")

ENlogReg_CN <- colnames(df_CBSA)
roc <- roc(as.formula(paste("Condition", "~", paste(ENlogReg_CN[-length(ENlogReg_CN)], collapse = "+"))), data = df_CBSA, family = 'binomial')
auc <- lapply(roc, function(i){as.numeric(i$auc)}) %>% data.frame() %>% t()
auc <- as.data.frame(auc)
colnames(auc) <- "AUC"
auc <- auc[order(-auc$AUC), , drop = FALSE]
auc

plot(roc$P01011, col = 'red')
legend(x = "bottomright", legend=c(paste("P01011", " (AUC = ", round(roc$P01011$auc,4), ")", sep = "")),
       col=c("red"), lty=1, cex=0.75, lwd = 2)

#####################################################
### CB vs SF Elastic Net Logistic Regression ########
df_CBSF <- df[df$Condition == "CB" | df$Condition == "SF",]
df_CBSF$Condition <- as.character(df_CBSF$Condition)
df_CBSF$Condition <- as.factor(df_CBSF$Condition)
df_CBSF$Condition <- relevel(df_CBSF$Condition, ref = "CB")
with(df_CBSF, table(Condition))

significantProteins.CBSF <- significantProteins[significantProteins$CB == "CB" | significantProteins$SF == "SF",]$Protein
significantProteins.CBSF
df_CBSF <- df_CBSF[c(significantProteins.CBSF, "Condition")]

ENlogReg_CBSF <- ENlogReg2(data = df_CBSF, k = 5)
ENlogReg_CBSF_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_CBSF, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_CBSF_Results, "ENlogReg_CBSF_Results.csv")

#####################################################
### SA vs SF Elastic Net Logistic Regression ########
df_SASF <- df[df$Condition == "SA" | df$Condition == "SF",]
df_SASF$Condition <- as.character(df_SASF$Condition)
df_SASF$Condition <- as.factor(df_SASF$Condition)
df_SASF$Condition <- relevel(df_SASF$Condition, ref = "SA")
with(df_SASF, table(Condition))

significantProteins.SASF <- significantProteins[significantProteins$SA == "SA" | significantProteins$SF == "SF",]$Protein
significantProteins.SASF
df_SASF <- df_SASF[c(significantProteins.SASF, "Condition")]

ENlogReg_SASF <- ENlogReg2(data = df_SASF, k = 5)
ENlogReg_SASF_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_SASF, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_SASF_Results, "ENlogReg_SASF_Results.csv")

#####################################################
df <- df.Imputed
df$Condition <- as.character(df$Condition)
df <- df[df$Condition == "HC" | df$Condition == "NSFM" | df$Condition == "CB" | df$Condition == "SA" | df$Condition == "SF",]
df$Condition <- gsub("SF", "SFM", df$Condition)
df$Condition <- gsub("CB", "SFM", df$Condition)
df$Condition <- gsub("SA", "SFM", df$Condition)
df$Condition <- gsub("NSFMM", "NSFM", df$Condition)
df$Condition <- as.factor(df$Condition)
df$Condition <- relevel(df$Condition, ref = "HC")
df$Condition
with(df, table(Condition))
# Condition
# HC NSFM  SFM 
# 18    8   24

#####################################################
### HC vs SFM Elastic Net Logistic Regression ########
df_HCSFM <- df[df$Condition == "HC" | df$Condition == "SFM",]
df_HCSFM$Condition <- as.character(df_HCSFM$Condition)
df_HCSFM$Condition <- as.factor(df_HCSFM$Condition)
df_HCSFM$Condition <- relevel(df_HCSFM$Condition, ref = "HC")
with(df_HCSFM, table(Condition))

significantProteins.HCSFM <- significantProteins[significantProteins$SF == "SF" | significantProteins$CB == "CB" | significantProteins$SA == "SA",]$Protein
significantProteins.HCSFM
df_HCSFM <- df_HCSFM[c(significantProteins.HCSFM, "Condition")]

ENlogReg_HCSFM <- ENlogReg(data = df_HCSFM, k = 10)
ENlogReg_HCSFM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCSFM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCSFM_Results, "ENlogReg_HCSFM_Results.csv")

#####################################################
### NSFM vs SFM Elastic Net Logistic Regression ########
df_NSFMSFM <- df[df$Condition == "NSFM" | df$Condition == "SFM",]
df_NSFMSFM$Condition <- as.character(df_NSFMSFM$Condition)
df_NSFMSFM$Condition <- as.factor(df_NSFMSFM$Condition)
df_NSFMSFM$Condition <- relevel(df_NSFMSFM$Condition, ref = "NSFM")
with(df_NSFMSFM, table(Condition))

significantProteins.NSFMSFM <- significantProteins[significantProteins$SF == "SF" | significantProteins$CB == "CB" | significantProteins$SA == "SA" | significantProteins$NSFM == "NSFM",]$Protein
significantProteins.NSFMSFM
df_NSFMSFM <- df_NSFMSFM[c(significantProteins.NSFMSFM, "Condition")]

ENlogReg_NSFMSFM <- ENlogReg2(data = df_NSFMSFM, k = 5)
ENlogReg_NSFMSFM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_NSFMSFM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_NSFMSFM_Results, "ENlogReg_NSFMSFM_Results.csv")

#########################################
#########################################
### Vivax Malaria - HC vs NSVM vs SVM ####
df <- df.Imputed
# Condition
# HC   CB  CCB  CSA   DF NSFM NSVM   SA   SD   SF  SVM 
# 18    6    4    4    4    8   12    6   12   12   20
df$Condition <- as.character(df$Condition)
df <- df[df$Condition == "HC" | df$Condition == "NSVM" | df$Condition == "SVM",]
df$Condition <- as.factor(df$Condition)
df$Condition <- relevel(df$Condition, ref = "HC")
df$Condition
with(df, table(Condition))
# Condition
# HC NSVM  SVM 
# 18   12   20  

#####################################################
### HC vs NSVM Elastic Net Logistic Regression ####
df_HCNSVM <- df[df$Condition != "SVM",]
df_HCNSVM$Condition <- as.character(df_HCNSVM$Condition)
df_HCNSVM$Condition <- as.factor(df_HCNSVM$Condition)
df_HCNSVM$Condition <- relevel(df_HCNSVM$Condition, ref = "HC")
with(df_HCNSVM, table(Condition))

significantProteins.HCNSVM <- significantProteins[significantProteins$VM == "VM",]$Protein
significantProteins.HCNSVM
df_HCNSVM <- df_HCNSVM[c(significantProteins.HCNSVM, "Condition")]

ENlogReg_HCNSVM <- ENlogReg(data = df_HCNSVM, k = 10)
ENlogReg_HCNSVM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCNSVM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCNSVM_Results, "ENlogReg_HCNSVM_Results.csv")

#####################################################
### HC vs SVM Elastic Net Logistic Regression ####
df_HCSVM <- df[df$Condition != "NSVM",]
df_HCSVM$Condition <- as.character(df_HCSVM$Condition)
df_HCSVM$Condition <- as.factor(df_HCSVM$Condition)
df_HCSVM$Condition <- relevel(df_HCSVM$Condition, ref = "HC")
with(df_HCSVM, table(Condition))

significantProteins.HCSVM <- significantProteins[significantProteins$VM == "VM",]$Protein
significantProteins.HCSVM
df_HCSVM <- df_HCSVM[c(significantProteins.HCSVM, "Condition")]

ENlogReg_HCSVM <- ENlogReg(data = df_HCSVM, k = 10)
ENlogReg_HCSVM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCSVM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCSVM_Results, "ENlogReg_HCSVM_Results.csv")

#####################################################
### NSVM vs SVM Elastic Net Logistic Regression ####
df_NSVMSVM <- df[df$Condition != "HC",]
df_NSVMSVM$Condition <- as.character(df_NSVMSVM$Condition)
df_NSVMSVM$Condition <- as.factor(df_NSVMSVM$Condition)
df_NSVMSVM$Condition <- relevel(df_NSVMSVM$Condition, ref = "NSVM")
with(df_NSVMSVM, table(Condition))

significantProteins.NSVMSVM <- significantProteins[significantProteins$VM == "VM",]$Protein
significantProteins.NSVMSVM
df_NSVMSVM <- df_NSVMSVM[c(significantProteins.NSVMSVM, "Condition")]

ENlogReg_NSVMSVM <- ENlogReg(data = df_NSVMSVM, k = 10)
ENlogReg_NSVMSVM_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_NSVMSVM, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_NSVMSVM_Results, "ENlogReg_NSVMSVM_Results.csv")

###################################################################
###################################################################
### Dengue Severity - HC vs DF vs SD ####
df <- df.Imputed
# Condition
# HC   CB  CCB  CSA   DF NSFM NSVM   SA   SD   SF  SVM 
# 18    6    4    4    4    8   12    6   12   12   20
df$Condition <- as.character(df$Condition)
df <- df[df$Condition == "HC" | df$Condition == "DF" | df$Condition == "SD",]
df$Condition <- as.factor(df$Condition)
df$Condition <- relevel(df$Condition, ref = "HC")
df$Condition
with(df, table(Condition))
# Condition
# HC DF SD 
# 18  4 12 

#######################################################
### HC vs DF Elastic Net Logistic Regression ########
df_HCDF <- df[df$Condition == "HC" | df$Condition == "DF",]
df_HCDF$Condition <- as.character(df_HCDF$Condition)
df_HCDF$Condition <- as.factor(df_HCDF$Condition)
df_HCDF$Condition <- relevel(df_HCDF$Condition, ref = "HC")
with(df_HCDF, table(Condition))

significantProteins.HCDF <- significantProteins[significantProteins$DF == "DF",]$Protein
significantProteins.HCDF
df_HCDF <- df_HCDF[c(significantProteins.HCDF, "Condition")]

ENlogReg_HCDF <- ENlogReg3(data = df_HCDF, k = 4)
ENlogReg_HCDF_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCDF, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCDF_Results, "ENlogReg_HCDF_Results.csv")

#######################################################
### HC vs SD Elastic Net Logistic Regression ########
df_HCSD <- df[df$Condition == "HC" | df$Condition == "SD",]
df_HCSD$Condition <- as.character(df_HCSD$Condition)
df_HCSD$Condition <- as.factor(df_HCSD$Condition)
df_HCSD$Condition <- relevel(df_HCSD$Condition, ref = "HC")
with(df_HCSD, table(Condition))

significantProteins.HCSD <- significantProteins[significantProteins$DF == "DF",]$Protein
significantProteins.HCSD
df_HCSD <- df_HCSD[c(significantProteins.HCSD, "Condition")]

ENlogReg_HCSD <- ENlogReg(data = df_HCSD, k = 10)
ENlogReg_HCSD_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_HCSD, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_HCSD_Results, "ENlogReg_HCSD_Results.csv")

#######################################################
### DF vs SD Elastic Net Logistic Regression ########
df_DFSD <- df[df$Condition == "DF" | df$Condition == "SD",]
df_DFSD$Condition <- as.character(df_DFSD$Condition)
df_DFSD$Condition <- as.factor(df_DFSD$Condition)
df_DFSD$Condition <- relevel(df_DFSD$Condition, ref = "DF")
with(df_DFSD, table(Condition))

significantProteins.DFSD <- significantProteins[significantProteins$DF == "DF",]$Protein
significantProteins.DFSD
df_DFSD <- df_DFSD[c(significantProteins.DFSD, "Condition")]

ENlogReg_DFSD <- ENlogReg3(data = df_DFSD, k = 4)
ENlogReg_DFSD_Results <- lapply(c('NumTrain','NumTest','TrainC','TestC','alpha','lambda','varImp','coef','kappa','cMbyClass'), function(i){ aggregateENModel(ENlogReg_DFSD, i) }) %>% do.call(rbind,.)
write.csv(ENlogReg_DFSD_Results, "ENlogReg_DFSD_Results.csv")
