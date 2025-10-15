classifyMATRIX <- function(features.sl, analysis = "", classes = NULL, class1 = "Healthy", class2 = "Cancer", testIND = NULL, method = "glmnet"){
  
  if(method=="gbm"){
    tuneGrid <- data.frame(n.trees=150,
                           interaction.depth=3,
                           shrinkage=0.1,
                           n.minobsinnode=10)
  }else if(method =="glmnet"){
    tuneGrid <- expand.grid(.alpha = 0.5, .lambda = seq(0, 0.05, by = 0.01)) 
  }else if(method =="svmLinear"){
    tuneGrid <- expand.grid(C = seq(0, 2, length = 20))
  }
  
  samples <- rownames(features.sl)
  features.sl$Class <- classes
  features.sl$Class <- factor(features.sl$Class, levels = c(unique(features.sl$Class)[unique(features.sl$Class)!=class1],unique(features.sl$Class)[unique(features.sl$Class)==class1]))
  
  if(!is.null(testIND)){
    samplesTRAIN <- samples[setdiff(1:length(samples),testIND)]
    samplesTEST <- samples[testIND]
    features.sl_TRAIN <- features.sl[samplesTRAIN,]
    features.sl_TEST <- features.sl[samplesTEST,]
    features.sl <- features.sl_TRAIN
  }else{
    features.sl <- features.sl[samples,]
  }
  
  summary.df <- data.frame(row.names = samples)
  summary.df$Class = classes
  
  if(!is.null(classes)) summary.df$Class <- classes
  
  features.sl$Class <- ifelse(features.sl$Class == class1, class1, class2)
  
  ctrl <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 100,
                       verboseIter = FALSE,
                       savePredictions=TRUE,
                       classProbs=TRUE,
                       summaryFunction = twoClassSummary,
                       allowParallel=TRUE)
  
  set.seed(1234)
  
  library(parallel) 
  no_cores <- 20
  
  library(doParallel)
  cl <- makePSOCKcluster(no_cores)
  registerDoParallel(cl)
  
  library(caret)
  
  model_gbm <- train(Class ~ .,
                     data = features.sl,
                     method = method,
                     tuneGrid=tuneGrid,
                     trControl = ctrl)
  
  stopCluster(cl)
  registerDoSEQ()
  
  if(!is.null(testIND)) {
    pred.tbl <- predict.train(model_gbm, newdata = features.sl_TEST, type = "prob")
    rownames(pred.tbl) <- rownames(features.sl_TEST)
    
    pred.tbl$Class <- features.sl_TEST$Class
    pred.tbl
  }else{
    
    if(sum(is.na(model_gbm$pred[,3]))>0){
      print("NA VALUE PRED")
      model_gbm$pred <- model_gbm$pred[!is.na(model_gbm$pred[,3]),]
    }
    
    pred.tbl <- as.data.frame(model_gbm$pred %>% group_by(rowIndex) %>% dplyr::summarize(Score=mean(.data[[class2]]), obs=obs[1]))
    rownames(pred.tbl) <- rownames(features.sl)
    pred.tbl$rowIndex <- NULL
    
  }
  
  pred.tbl
  
}