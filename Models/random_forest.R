library(ggplot2)
library(data.table)
library(caret)
library(scales)
library(randomForest)
library(SDMTools)

set.seed(0)

#----------- Read data, split out validation set and split out labels
train <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')
test <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/test.csv')

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


#-------- Fit Random Forest model
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

g <- ggplot(data=importance, aes(x=Feature, y=MeanDecreaseGini, fill=Feature)) + geom_bar(stat='identity')
g <- g + labs(x='Feature', y='Mean Decrease Gini Index') + coord_flip() + guides(fill=FALSE)
g + scale_x_discrete(limits = rev(levels(importance$Feature)))

# # Save importance features for later use
# important_features <- as.character(importance$Feature)
# important_features <- c(important_features, 'name', 'sequence', 'class')
# save(important_features, file = "~/Desktop/CSML/bioinformatics/coursework/important_features.RData")


# #------ Assess performance on test set
# X_train <- rbind(X_train, X_val)
# y_train <- factor(c(as.character(y_train), as.character(y_val)), level=c('cyto', 'mito', 'nuclear', 'secreted'))
# train_set <- X_train
# train_set$y <- y_train
# rf_model <- randomForest(y ~ ., data=train_set, ntree = 100, mtry = 9, nodesize=1)
# 
# # Accuracy
# test_pred <- predict(rf_model, newdata=X_test)
# print("Test Accuracy:")
# print(sum(test_pred==y_test)/length(y_test))
# 
# # Confusion matrix
# table(y_test, test_pred)

