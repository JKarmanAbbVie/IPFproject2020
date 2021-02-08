######################Figure 8.##########################
library(ISLR)
library(caret)
library(MLeval)
options(stringsAsFactors = F)

#celltype-based classifiers
gse47460_celltypegsva_consensus <- merge(consensusdf_gse47460, t(gse47460_celltypegsva_table), by.x = 0, by.y = 0)
rownames(gse47460_celltypegsva_consensus) <- gse47460_celltypegsva_consensus$Row.names
gse47460_celltypegsva_consensus$Row.names <- NULL
gse47460_celltypegsva_consensus$Geo <- NULL
gse47460_celltypegsva_consensus$consensusclass <- paste0('cluster_', gse47460_celltypegsva_consensus$consensusclass)
gse47460_celltypegsva_consensus$consensusclass <- as.factor(gse47460_celltypegsva_consensus$consensusclass)
set.seed(1)
indxTrain <- createDataPartition(y = gse47460_celltypegsva_consensus$consensusclass,p = 0.7,list = FALSE)
training <- gse47460_celltypegsva_consensus[indxTrain,]
testing <- gse47460_celltypegsva_consensus[-indxTrain,]
set.seed(400)
ctrl <- trainControl(method="repeatedcv", repeats = 5 , classProbs=TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
svmFit <- train(consensusclass ~ ., data = training, method = "svmLinear", trControl = ctrl, preProcess = c("center","scale"))
svmPredict <- predict(svmFit,newdata = testing )
library(pROC)
svmPredict <- predict(svmFit,newdata = testing , type="prob")
svmROC <- roc(as.factor(testing$consensusclass),svmPredict[,"cluster_1"])
svmimportance <- varImp(svmFit)
library(gbm)
gbmFit <- train(consensusclass ~ ., data = training, method = "gbm", trControl = ctrl, preProcess = c("center","scale"))
gbmPredict <- predict(gbmFit,newdata = testing , type="prob")
gbmROC <- roc(as.factor(testing$consensusclass),gbmPredict[,"cluster_1"])
gbmimportance <- varImp(gbmFit, scale = F)
glmnetFit <- train(consensusclass ~ ., data = training, method = "glmnet", trControl = ctrl, preProcess = c("center","scale"))
glmnetPredict <- predict(glmnetFit,newdata = testing , type="prob")
glmnetROC <- roc(as.factor(testing$consensusclass),glmnetPredict[,"cluster_1"])
glmnetimportance <- varImp(glmnetFit, scale = F)
#Figure 8B:
print(glmnetimportance)

library(mlbench)
res <- evalm(list(svmFit,gbmFit,glmnetFit),gnames=c('svm','gbm','glmnet'))

#gene-based classifiers
gse47460_phenom_consensus_ipf <- as.data.frame(merge(consensusdf_gse47460, gse47460_phenom, by.x = 'Geo', by.y = 0))
gse47460_phenom_consensus_ipf$Geo <- NULL
gse47460_phenom_consensus_ipf$consensusclass <- paste0('cluster_', gse47460_phenom_consensus_ipf$consensusclass)
gse47460_phenom_consensus_ipf$consensusclass <- as.factor(gse47460_phenom_consensus_ipf$consensusclass)
gse47460_phenom_consensus_ipf <- gse47460_phenom_consensus_ipf[,c(1, 15:19010)]
mads_gse47460=apply(gse47460_phenom_consensus_ipf[,2:18997],2,mad)
gse474760_mostvariable <- names(sort(mads_gse47460, decreasing = T)[1:10000])
gse47460_phenom_consensus_ipf <- gse47460_phenom_consensus_ipf[,c('consensusclass', gse474760_mostvariable)]

set.seed(1)
indxTrain <- createDataPartition(y = gse47460_celltypegsva_consensus$consensusclass,p = 0.7,list = FALSE)
training <- gse47460_celltypegsva_consensus[indxTrain,]
testing <- gse47460_celltypegsva_consensus[-indxTrain,]
set.seed(400)
ctrl <- trainControl(method="repeatedcv", repeats = 5 , classProbs=TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)
svmFit <- train(consensusclass ~ ., data = training, method = "svmLinear", trControl = ctrl, preProcess = c("center","scale"))
svmPredict <- predict(svmFit,newdata = testing )
library(pROC)
svmPredict <- predict(svmFit,newdata = testing , type="prob")
svmROC <- roc(as.factor(testing$consensusclass),svmPredict[,"cluster_1"])
svmimportance <- varImp(svmFit)
library(gbm)
gbmFit <- train(consensusclass ~ ., data = training, method = "gbm", trControl = ctrl, preProcess = c("center","scale"))
gbmPredict <- predict(gbmFit,newdata = testing , type="prob")
gbmROC <- roc(as.factor(testing$consensusclass),gbmPredict[,"cluster_1"])
gbmimportance <- varImp(gbmFit, scale = F)
glmnetFit <- train(consensusclass ~ ., data = training, method = "glmnet", trControl = ctrl, preProcess = c("center","scale"))
glmnetPredict <- predict(glmnetFit,newdata = testing , type="prob")
glmnetROC <- roc(as.factor(testing$consensusclass),glmnetPredict[,"cluster_1"])
glmnetimportance <- varImp(glmnetFit, scale = F)
#Figure 8C:
res <- evalm(list(svmFit,gbmFit,glmnetFit),gnames=c('svm','gbm','glmnet'))

#Recursive feature elimination
control <- rfeControl(functions=caretFuncs, method="cv", number=5)
results_rfe <- rfe(gse47460_phenom_consensus_ipf[,2:10001], gse47460_celltypegsva_consensus[,1], sizes=c(1:50), rfeControl=control)
print(results_rfe)
predictors(results_rfe)
plot(results_rfe, type=c("g", "o"))

#Figure 8D: same code below repeated for each gene, example of FOXJ1 shown
ggplot(gse47460_phenom_consensus, aes(consensusclass, FOXJ1, fill = consensusclass)) + 
  geom_boxplot(width=0.5) + geom_jitter(color="black", shape=16, position=position_jitter(0.2)) + theme_bw() + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

dunnTest(FOX11~consensusclass, data = gse47460_phenom_consensus)
#end
