rm(list = ls())

require(readr)
require(ggplot2)
require(dplyr)
require(tidyr)
require(caret)
require(corrplot)
require(Hmisc)
require(parallel)
require(doParallel)
require(ggthemes)
library(ggsci)


set.seed(52)

load("data/step6_all_data.Rdata")

cg = c("TTC3", "MGST1", 'ACTR3')

data = all_data[-(1:5)][colnames(all_data[-(1:5)]) %in% cg]
#data = log2(data + 1)
data$group = all_data$group
data$group = ifelse(data$group == "NCG_low", "LRG", "HRG")
data$group = factor(data$group, levels = c("LRG", "HRG"))

### separate dataset into training and testing sets
sample_Index <- createDataPartition(data$group, p=0.7,list=FALSE)
voice_Train <- data[sample_Index,]
voice_Test <- data[-sample_Index,]

### preprocess factors for further modeling
pp <- preProcess(voice_Train,method=c("scale","center","pca"))
voice_Train <- predict(pp,voice_Train)
voice_Test <- predict(pp,voice_Test)

### define formula
model_Formula <- group~PC1+PC2+PC3

###set cross-validation parameters
modelControl <- trainControl(method="repeatedcv",number=5,
                             repeats=5,allowParallel=TRUE, 
                             classProbs = TRUE)


### logistic regression
glm_Model <- train(model_Formula,
                   data=voice_Train,
                   method="glm",
                   trControl=modelControl)

#glm_Coefficients <- coef(glm_Model$finalModel)

### linear discrimant analysis
lda_Model <- train(model_Formula,
                   data=voice_Train,
                   method="lda",
                   trControl=modelControl)

### random forrest
rf_Model <- train(model_Formula,
                  data=voice_Train,
                  method="rf",
                  trControl=modelControl,
                  ntrees=500)

### Naive Bayes
nb_Model <- train(model_Formula,
                  data=voice_Train,
                  method="nb",
                  trControl=modelControl)

### SVM
svm_Model <- train(model_Formula,
                  data=voice_Train,
                  method="svmLinear",
                  trControl=modelControl,
                  )

### gbm
gbm_Model <- train(model_Formula,
                   data=voice_Train,
                   method="gbm",
                   trControl=modelControl,
)

### xgboost
xgboost_Model <- train(model_Formula,
                   data=voice_Train,
                   method="xgbLinear",
                   trControl=modelControl,
)

### nnet
nnet_Model <- train(model_Formula,
                       data=voice_Train,
                       method="nnet",
                       trControl=modelControl,
)

### knn
knn_Model <- train(model_Formula,
                    data=voice_Train,
                    method="knn",
                    trControl=modelControl,
)

### rpart
rpart_Model <- train(model_Formula,
                   data=voice_Train,
                   method="C5.0",
                   trControl=modelControl,
)

### AdaBoost
AdaBoost_Model <- train(model_Formula,
                     data=voice_Train,
                     method="AdaBoost.M1",
                     trControl=modelControl,
)

library(pROC)

a = voice_Test
a$glm_Prob <- predict(glm_Model, voice_Test, type = "prob")[, 2]
a$lda_Prob <- predict(lda_Model, voice_Test, type = "prob")[, 2]
a$rf_Prob <- predict(rf_Model, voice_Test, type = "prob")[, 2]
a$nb_Prob <- predict(nb_Model, voice_Test, type = "prob")[, 2]
a$svm_Prob <- predict(svm_Model, voice_Test, type = "prob")[, 2]
a$gbm_Prob <- predict(gbm_Model, voice_Test, type = "prob")[, 2]
a$xgboost_Prob <- predict(xgboost_Model, voice_Test, type = "prob")[, 2]
a$nnet_Prob <- predict(nnet_Model, voice_Test, type = "prob")[, 2]
a$knn_Prob <- predict(knn_Model, voice_Test, type = "prob")[, 2]
a$rpart_Prob <- predict(rpart_Model, voice_Test, type = "prob")[, 2]
a$AdaBoost_Prob <- predict(AdaBoost_Model, voice_Test, type = "prob")[, 2]

save(a, voice_Test, voice_Train, file = "./data/CESC_Predict_modle.Rdata")
load("./data/CESC_Predict_modle.Rdata")

colors=pal_futurama()(11)

pdf(file = "./Figure/step6_5.pdf", onefile = F, width = 8)
plot.roc(
  a$group,              
  a$glm_Prob,       
  col = colors[1],            
  percent = TRUE,             
  lwd = 2,                   
  print.auc = TRUE,           
  print.auc.cex = 1,          
  print.auc.pattern = "Logistic Regression: AUC=%.1f%%",      
  print.auc.y = 50,
  main = "TCGA-CESC"
)

plot.roc(a$group,    
         a$lda_Prob,
         col = colors[2],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Linear Regression: AUC=%.1f%%",
         print.auc.y = 45,
         add = T)

plot.roc(a$group,    
         a$rf_Prob,
         col = colors[3],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Random Forest: AUC=%.1f%%",
         print.auc.y = 40,
         add = T)

plot.roc(a$group,    
         a$nb_Prob,
         col = colors[4],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Naive Bayes: AUC=%.1f%%",
         print.auc.y = 35,
         add = T)

plot.roc(a$group,    
         a$svm_Prob,
         col = colors[5],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Support Vector Machine: AUC=%.1f%%",
         print.auc.y = 30,
         add = T)

plot.roc(a$group,    
         a$gbm_Prob,
         col = colors[6],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Gradient Boosting Machine: AUC=%.1f%%",
         print.auc.y = 25,
         add = T)

plot.roc(a$group,    
         a$xgboost_Prob,
         col = colors[7],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Extreme Gradient Boosting: AUC=%.1f%%",
         print.auc.y = 20,
         add = T)

plot.roc(a$group,    
         a$nnet_Prob,
         col = colors[8],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Artificial Neural Network: AUC=%.1f%%",
         print.auc.y = 15,
         add = T)

plot.roc(a$group,    
         a$knn_Prob,
         col = colors[9],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "k-Nearest Neighbors: AUC=%.1f%%",
         print.auc.y = 10,
         add = T)

plot.roc(a$group,    
         a$rpart_Prob,
         col = colors[10],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Decision Tree: AUC=%.1f%%",
         print.auc.y = 5,
         add = T)

plot.roc(a$group,    
         a$AdaBoost_Prob,
         col = colors[11],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Adaptive Boosting: AUC=%.1f%%",
         print.auc.y = 0,
         add = T)
dev.off()

### GSE56363  != hg38
#rm(list = ls())

cg = c("TTC3", "MGST1", 'ACTR3')

GEE1359941 = read.table("./data/GSE56363_series_matrix.txt.gz", fill = T, header = T, comment.char = "!")
GEE1359941_anno = read.table("./data/GPL4133_old_annotations.txt.gz", fill = T, comment.char = "#", skip = 1, header = T)

GEE1359941_anno = GEE1359941_anno[c(1,9:14)]
colnames(GEE1359941_anno)[1] = "ID_REF"
GEE1359941$ID_REF = as.character(GEE1359941$ID_REF)

GEE1359941_sub = GEE1359941_anno %>% 
  inner_join(GEE1359941, by="ID_REF") %>% 
  filter(GENE_SYMBOL %in% cg)

GEE1359941_sub = GEE1359941_sub[-c(1:2, 4:7)]
nor_data = GEE1359941_sub

nor_data$mean_expression <- rowMeans(nor_data[2:22])
nor_data <- nor_data[order(nor_data$GENE_SYMBOL, -nor_data$mean_expression), ]
nor_data <- nor_data[!duplicated(nor_data$GENE_SYMBOL), ]
rownames(nor_data) = nor_data$GENE_SYMBOL
nor_data = nor_data[2:22]

nor_data = as.data.frame(t(nor_data))

nor_data$group = c(rep("LRG", 12), rep("HRG", 9))
nor_data$group = factor(nor_data$group, levels = c("LRG", "HRG"))
voice_Test <- nor_data
voice_Test <- predict(pp,voice_Test)

a = voice_Test
a$glm_Prob <- predict(glm_Model, voice_Test, type = "prob")[, 1]
a$lda_Prob <- predict(lda_Model, voice_Test, type = "prob")[, 1]
a$rf_Prob <- predict(rf_Model, voice_Test, type = "prob")[, 1]
a$nb_Prob <- predict(nb_Model, voice_Test, type = "prob")[, 1]
a$svm_Prob <- predict(svm_Model, voice_Test, type = "prob")[, 1]
a$gbm_Prob <- predict(gbm_Model, voice_Test, type = "prob")[, 1]
a$xgboost_Prob <- predict(xgboost_Model, voice_Test, type = "prob")[, 1]
a$nnet_Prob <- predict(nnet_Model, voice_Test, type = "prob")[, 1]
a$knn_Prob <- predict(knn_Model, voice_Test, type = "prob")[, 1]
a$rpart_Prob <- predict(rpart_Model, voice_Test, type = "prob")[, 1]
a$AdaBoost_Prob <- predict(AdaBoost_Model, voice_Test, type = "prob")[, 1]

pdf(file = "./Figure/step6_6.pdf", onefile = F, width = 8)
plot.roc(
  a$group,              
  a$glm_Prob,       
  col = colors[1],            
  percent = TRUE,             
  lwd = 2,                   
  print.auc = TRUE,           
  print.auc.cex = 1,          
  print.auc.pattern = "Logistic Regression: AUC=%.1f%%",      
  print.auc.y = 50,
  main = "GSE56363"
)

plot.roc(a$group,    
         a$lda_Prob,
         col = colors[2],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Linear Regression: AUC=%.1f%%",
         print.auc.y = 45,
         add = T)

plot.roc(a$group,    
         a$rf_Prob,
         col = colors[3],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Random Forest: AUC=%.1f%%",
         print.auc.y = 40,
         add = T)

plot.roc(a$group,    
         a$nb_Prob,
         col = colors[4],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Naive Bayes: AUC=%.1f%%",
         print.auc.y = 35,
         add = T)

plot.roc(a$group,    
         a$svm_Prob,
         col = colors[5],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Support Vector Machine: AUC=%.1f%%",
         print.auc.y = 30,
         add = T)

plot.roc(a$group,    
         a$gbm_Prob,
         col = colors[6],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Gradient Boosting Machine: AUC=%.1f%%",
         print.auc.y = 25,
         add = T)

plot.roc(a$group,    
         a$xgboost_Prob,
         col = colors[7],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Extreme Gradient Boosting: AUC=%.1f%%",
         print.auc.y = 20,
         add = T)

plot.roc(a$group,    
         a$nnet_Prob,
         col = colors[8],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Artificial Neural Network: AUC=%.1f%%",
         print.auc.y = 15,
         add = T)

plot.roc(a$group,    
         a$knn_Prob,
         col = colors[9],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "k-Nearest Neighbors: AUC=%.1f%%",
         print.auc.y = 10,
         add = T)

plot.roc(a$group,    
         a$rpart_Prob,
         col = colors[10],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Decision Tree: AUC=%.1f%%",
         print.auc.y = 5,
         add = T)

plot.roc(a$group,    
         a$AdaBoost_Prob,
         col = colors[11],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Adaptive Boosting: AUC=%.1f%%",
         print.auc.y = 0,
         add = T)
dev.off()

### GSE168009
GSE168 = read.table("./data/GSE168009_Raw_count.txt.gz", fill = T, header = T)
GSE168 = as.data.frame(t(GSE168))

GSE168 = GSE168[colnames(GSE168) %in% cg]
GSE168 = GSE168[-c(6,8),]
GSE168$group = c(rep("HRG", 4), rep("LRG", 3))
GSE168$group = factor(GSE168$group, levels = c("LRG", "HRG"))

voice_Test <- GSE168
voice_Test <- predict(pp,voice_Test)

a = voice_Test
a$glm_Prob <- predict(glm_Model, voice_Test, type = "prob")[, 2]
a$lda_Prob <- predict(lda_Model, voice_Test, type = "prob")[, 2]
a$rf_Prob <- predict(rf_Model, voice_Test, type = "prob")[, 2]
a$nb_Prob <- predict(nb_Model, voice_Test, type = "prob")[, 2]
a$svm_Prob <- predict(svm_Model, voice_Test, type = "prob")[, 2]
a$gbm_Prob <- predict(gbm_Model, voice_Test, type = "prob")[, 2]
a$xgboost_Prob <- predict(xgboost_Model, voice_Test, type = "prob")[, 2]
a$nnet_Prob <- predict(nnet_Model, voice_Test, type = "prob")[, 2]
a$knn_Prob <- predict(knn_Model, voice_Test, type = "prob")[, 2]
a$rpart_Prob <- predict(rpart_Model, voice_Test, type = "prob")[, 2]
a$AdaBoost_Prob <- predict(AdaBoost_Model, voice_Test, type = "prob")[, 2]

pdf(file = "./Figure/step6_7.pdf", onefile = F, width = 8)
plot.roc(
  a$group,              
  a$glm_Prob,       
  col = colors[1],            
  percent = TRUE,             
  lwd = 2,                   
  print.auc = TRUE,           
  print.auc.cex = 1,          
  print.auc.pattern = "Logistic Regression: AUC=%.1f%%",      
  print.auc.y = 50,
  main = "GSE168009"
)

plot.roc(a$group,    
         a$lda_Prob,
         col = colors[2],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Linear Regression: AUC=%.1f%%",
         print.auc.y = 45,
         add = T)

plot.roc(a$group,    
         a$rf_Prob,
         col = colors[3],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Random Forest: AUC=%.1f%%",
         print.auc.y = 40,
         add = T)

plot.roc(a$group,    
         a$nb_Prob,
         col = colors[4],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Naive Bayes: AUC=%.1f%%",
         print.auc.y = 35,
         add = T)

plot.roc(a$group,    
         a$svm_Prob,
         col = colors[5],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Support Vector Machine: AUC=%.1f%%",
         print.auc.y = 30,
         add = T)

plot.roc(a$group,    
         a$gbm_Prob,
         col = colors[6],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Gradient Boosting Machine: AUC=%.1f%%",
         print.auc.y = 25,
         add = T)

plot.roc(a$group,    
         a$xgboost_Prob,
         col = colors[7],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Extreme Gradient Boosting: AUC=%.1f%%",
         print.auc.y = 20,
         add = T)

plot.roc(a$group,    
         a$nnet_Prob,
         col = colors[8],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Artificial Neural Network: AUC=%.1f%%",
         print.auc.y = 15,
         add = T)

plot.roc(a$group,    
         a$knn_Prob,
         col = colors[9],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "k-Nearest Neighbors: AUC=%.1f%%",
         print.auc.y = 10,
         add = T)

plot.roc(a$group,    
         a$rpart_Prob,
         col = colors[10],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Decision Tree: AUC=%.1f%%",
         print.auc.y = 5,
         add = T)

plot.roc(a$group,    
         a$AdaBoost_Prob,
         col = colors[11],
         percent = T,
         lwd = 2,
         print.auc = T,
         print.auc.cex = 1,
         print.auc.pattern = "Adaptive Boosting: AUC=%.1f%%",
         print.auc.y = 0,
         add = T)
dev.off()

