library(ggplot2)
library(data.table)
library(caret)
library(scales)
library(randomForest)
library(SDMTools)
library(pROC)
library(klaR)
library(xtable)

set.seed(0)

#----------- Part 1: Read data, split out validation set and split out labels
train <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')
test <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/test.csv')
blind_test <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/blind_test.csv')

# # Use important features only
# load(file = "~/Desktop/CSML/bioinformatics/coursework/important_features.RData")
# train <- train[, important_features]
# test <- test[, important_features]

smp_size <- floor(0.75 * nrow(train))
train_ind <- sample(seq_len(nrow(train)), size = smp_size)
train_set <- train[train_ind, ]
val_set <- train[-train_ind, ]

y_train <- factor(train_set$class, level=c('cyto', 'mito', 'nuclear', 'secreted'))
names_train <- train_set$name
X_train <- subset(train_set, select=-c(name, class, sequence))

y_val <- factor(val_set$class, level=c('cyto', 'mito', 'nuclear', 'secreted'))
names_val <- val_set$name
X_val <- subset(val_set, select=-c(name, class, sequence))

y_test <- factor(test$class, level=c('cyto', 'mito', 'nuclear', 'secreted'))
names_test <- test$name
X_test <- subset(test, select=-c(name, class, sequence))


#-------- Part 2: Fit Random Forest model
train_set <- X_train
train_set$y <- y_train
rf_model <- randomForest(y ~ ., data=train_set, ntree = 500, mtry = 9, nodesize=1)

# Train and validation set predictions
train_pred <- predict(rf_model)
print("Training Accuracy:")
print(sum(train_pred==y_train)/length(y_train))
table(y_train, train_pred)

val_pred <- predict(rf_model, newdata=X_val)
print("Validation Accuracy:")
print(sum(val_pred==y_val)/length(y_val))

# Feature importance
importance <- data.frame(rf_model$importance)
importance$Feature <- rownames(importance)
importance <- data.table(importance)[order(-MeanDecreaseGini), ]
setcolorder(importance, c('Feature', 'MeanDecreaseGini'))
importance <- importance[MeanDecreaseGini>=29]
importance$Feature <- factor(importance$Feature)
importance$Feature <- reorder(importance$Feature, -importance$MeanDecreaseGini)

# Plot and save as PNG file
g <- ggplot(data=importance, aes(x=Feature, y=MeanDecreaseGini, fill=Feature)) + geom_bar(stat='identity')
g <- g + labs(x='Feature', y='Mean Decrease Gini Index') + coord_flip() + guides(fill=FALSE)
g + scale_x_discrete(limits = rev(levels(importance$Feature)))

png(filename="~/Desktop/CSML/bioinformatics/coursework/figures/feature_importance.png")
g + scale_x_discrete(limits = rev(levels(importance$Feature)))
dev.off()

# # Save importance features for later use
# important_features <- as.character(importance$Feature)
# important_features <- c(important_features, 'name', 'sequence', 'class')
# save(important_features, file = "~/Desktop/CSML/bioinformatics/coursework/important_features.RData")


#------ Part 3: Assess performance on test set
X_train <- rbind(X_train, X_val)
y_train <- factor(c(as.character(y_train), as.character(y_val)), level=c('cyto', 'mito', 'nuclear', 'secreted'))
train_set <- X_train
train_set$y <- y_train
rf_model <- randomForest(y ~ ., data=train_set, ntree = 100, mtry = 9, nodesize=1)

# Accuracy
test_pred <- predict(rf_model, newdata=X_test)
print("Test Accuracy:")
print(sum(test_pred==y_test)/length(y_test))

# Confusion matrix
table(y_test, test_pred)
xtable(data.frame(table(y_test, test_pred))) # latex table format

# Plot ROC curves and compute AUC
test_pred_prob <- data.table(data.frame(predict(rf_model, newdata=X_test, type='prob')))
test_pred_prob$prediction <- test_pred
test_pred_prob$true_label <- y_test

class_prob_all <- data.frame()
auc_scores <- c()

for (class_name in c('cyto', 'mito', 'nuclear', 'secreted')) {
  class_prob <- subset(test_pred_prob, select=c(eval(as.symbol(class_name)), prediction, true_label))[order(-eval(as.symbol(class_name)))]
  
  # AUC score
  prob_scores <- as.numeric(unlist(subset(class_prob, select=eval(as.symbol(class_name)))))
  area <- auc(as.numeric(class_prob$true_label==class_name), prob_scores)[1]
  auc_scores <- c(auc_scores, area)
  
  # ROC curve
  class_prob$FP <- 0; class_prob$TP <- 0
  class_prob[class_prob$prediction == class_name & class_prob$true_label == class_name, ]$TP <- 1
  class_prob[class_prob$prediction == class_name & class_prob$true_label != class_name, ]$FP <- 1
  class_prob$pos <- 0; class_prob$neg <- 0
  class_prob[class_prob$true_label == class_name, ]$pos <- 1
  class_prob[class_prob$true_label != class_name, ]$neg <- 1
  class_prob$TPR <- cumsum(class_prob$TP) / sum(class_prob$TP)
  class_prob$FPR <- cumsum(class_prob$FP) /sum(class_prob$FP)
  class_prob <- subset(class_prob, select=c(TPR, FPR))
  class_prob$class <- class_name
  class_prob_all <- rbind(class_prob_all, class_prob)
}

class_prob_all$class <- factor(class_prob_all$class)
g <- ggplot(data=class_prob_all, aes(x=FPR, y=TPR, color=class)) + geom_line(size=0.8)
g <- g + labs(x='False positive rate', y='True positive rate', color='Protein location class')
g + theme(legend.position="top")

png(filename="~/Desktop/CSML/bioinformatics/coursework/figures/ROC_curve.png")
g + theme(legend.position="top")
dev.off()

# AUC table and mean AUC
print("Mean AUC:")
mean(auc_scores)
auc_scores <- data.frame(Class=c('cyto', 'mito', 'nuclear', 'secreted'), AUC_score=auc_scores)
xtable(auc_scores)  # latex output


##----- Part 4: Predictions on blind test set
protein_names <- blind_test$name
blind_test <- subset(blind_test, select=-c(name, sequence))
blind_test_pred <- data.table(data.frame(predict(rf_model, newdata=blind_test, type='prob')))
blind_test_pred$name <- protein_names
blind_test_pred <- melt(blind_test_pred, id.vars='name')
blind_test_pred <- blind_test_pred[, .SD[value==max(value)], by=name]
blind_test_pred <- blind_test_pred[order(name)]
names(blind_test_pred) <- c('Sequence Identifier', 'Class Prediction', 'Value')
xtable(blind_test_pred)