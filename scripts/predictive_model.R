normalizeMatrix <- function(feat_mtx, Cohort = ""){
  NUM_FEATURES <- ncol(feat_mtx)
  feat_mtx$Cohort = Cohort
  feat_mtx$Sample <- rownames(feat_mtx)
  feat_mtx_norm <- as.data.frame(feat_mtx %>% dplyr::group_by(Cohort) %>% dplyr::mutate_at(colnames(feat_mtx)[1:NUM_FEATURES], function(x) minMaxNorm(x, min(x), max(x))))
  feat_mtx_norm <- feat_mtx_norm[1:NUM_FEATURES]
  rownames(feat_mtx_norm) <- rownames(feat_mtx)
  feat_mtx_norm
}

classifyMATRIX <- function(feat_mtx, classes = NULL, Cohort = "", class1 = "Healthy", class2 = "Cancer", testIND = NULL, method = "glmnet", NUM_THREADS = 20){
  
  feat_mtx <- normalizeMatrix(feat_mtx, Cohort)
  
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
  
  samples <- rownames(feat_mtx)
  feat_mtx$Class <- classes
  feat_mtx$Class <- factor(feat_mtx$Class, levels = c(unique(feat_mtx$Class)[unique(feat_mtx$Class)!=class1],unique(feat_mtx$Class)[unique(feat_mtx$Class)==class1]))
  
  if(!is.null(testIND)){
    samplesTRAIN <- samples[setdiff(1:length(samples),testIND)]
    samplesTEST <- samples[testIND]
    feat_mtx_TRAIN <- feat_mtx[samplesTRAIN,]
    feat_mtx_TEST <- feat_mtx[samplesTEST,]
    feat_mtx <- feat_mtx_TRAIN
  }else{
    feat_mtx <- feat_mtx[samples,]
  }
  
  summary.df <- data.frame(row.names = samples)
  summary.df$Class = classes
  
  if(!is.null(classes)) summary.df$Class <- classes
  
  feat_mtx$Class <- ifelse(feat_mtx$Class == class1, class1, class2)
  
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
  no_cores <- NUM_THREADS
  
  library(doParallel)
  cl <- makePSOCKcluster(no_cores)
  registerDoParallel(cl)
  
  library(caret)
  
  model_gbm <- train(Class ~ .,
                     data = feat_mtx,
                     method = method,
                     tuneGrid=tuneGrid,
                     trControl = ctrl)
  
  stopCluster(cl)
  registerDoSEQ()
  
  if(!is.null(testIND)) {
    pred.tbl <- predict.train(model_gbm, newdata = feat_mtx_TEST, type = "prob")
    rownames(pred.tbl) <- rownames(feat_mtx_TEST)
    
    pred.tbl$Class <- feat_mtx_TEST$Class
    pred.tbl
  }else{
    
    if(sum(is.na(model_gbm$pred[,3]))>0){
      print("NA VALUE PRED")
      model_gbm$pred <- model_gbm$pred[!is.na(model_gbm$pred[,3]),]
    }
    
    pred.tbl <- as.data.frame(model_gbm$pred %>% group_by(rowIndex) %>% dplyr::summarize(Score=mean(.data[[class2]]), obs=obs[1]))
    rownames(pred.tbl) <- rownames(feat_mtx)
    pred.tbl$rowIndex <- NULL
    
  }
  
  pred.tbl
  
}