library(e1071)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(xgboost)
library(caret)
library(scales)

set.seed(0)

# Read and split the data
train <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')
test <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/test.csv')

# train <- train[, c("L_first50_count", "L_last50_count",  "A_first50_count", "A_last50_count",  "V_first50_count", 
#                    "V_last50_count", "molecular_weight", "sequence_length", "K_first50_count", "K_last50_count",  
#                    "F_first50_count", "F_last50_count", "second_helix","aromaticity", "E_count","isolectric_point",
#                    "K_count","V_count","E_first50_count", "E_last50_count",  "R_first50_count", "R_last50_count",
#                    "M_first50_count", "M_last50_count", "W_first50_count", "W_last50_count",  "F_count","S_count",
#                    "M_count","S_first50_count", 'name', 'sequence', 'class')]

smp_size <- floor(0.75 * nrow(train))
train_ind <- sample(seq_len(nrow(train)), size = smp_size)
train_set <- train[train_ind, ]
val_set <- train[-train_ind, ]

y_train <- factor(train_set$class, levels=c('cyto', 'mito', 'nuclear', 'secreted'))
names_train <- train_set$name
X_train <- subset(train_set, select=-c(name, class, sequence))

y_val <- factor(val_set$class, levels=c('cyto', 'mito', 'nuclear', 'secreted'))
names_val <- val_set$name
X_val <- subset(val_set, select=-c(name, class, sequence))

y_test <- factor(test$class, levels=c('cyto', 'mito', 'nuclear', 'secreted'))
names_test <- test$name
X_test <- subset(test, select=-c(name, class, sequence))



#-------------- Fit SVM model

## Function to calculate the accuracy for a given model and test dataset
label_list <- c('cyto', 'mito', 'nuclear', 'secreted')

predSVM <- function(model) {
  pred_prob <- predict(model, X_val, decision.values = TRUE, probability = TRUE)
  pred_df <- data.table(as.data.frame(attr(pred_prob, "probabilities")))
  label_list <- c('cyto', 'mito', 'nuclear', 'secreted')
  setcolorder(pred_df, label_list)
  y_hat_idx <- apply(pred_df, 1, function(x) max(which(x == max(x, na.rm = TRUE))))
  y_hat <- sapply(y_hat_idx, function(x) label_list[x])
  y_hat <- factor(y_hat, levels=label_list)
  
  # Calculate the prediction error
  results <- data.frame(actual=y_val, predicted=y_hat)
  results$mistake <- ifelse(results$actual!=results$predicted, 1, 0)
  return(1 - sum(results$mistake)/nrow(results))
}


# ## Use validation to choose a parameter C
# gamma_vals <- c(0.004, 0.005, 0.006)
# val_results <- data.table(cost_parameter=gamma_vals)
# val_results$val_error <- 9999
# 
# for (g in gamma_vals) {
#   print(paste('Running for gamma =', g))
#   model <- svm(X_train, y_train, cost=2, gamma=g, probability = TRUE)
#   error <- predSVM(model)
#   val_results[cost_parameter==g]$val_error <- error
# }
# val_results[order(-val_error)][1]$cost_parameter


model <- svm(X_train, y_train, cost=2, gamma=0.005, probability = TRUE)
predSVM(model)


# Feature importance
W <- subset(data.table(melt(t(model$coefs) %*% model$SV)), select=-Var1)
W$value <- abs(W$value)
W <- W[, lapply(.SD, mean), by=Var2]
# W$value <- abs(W$value)
W <- W[order(-value), ][1:30]
setnames(W, names(W), c('Feature', 'mean_coef_weight'))
W$Feature <- factor(W$Feature)
W$Feature <- reorder(W$Feature, -W$mean_coef_weight)

g <- ggplot(data=W, aes(x=Feature, y=mean_coef_weight, fill=Feature)) + geom_bar(stat='identity')
g <- g + labs(x='Feature', y='Mean Coef Weight') + coord_flip() + guides(fill=FALSE)
g + scale_x_discrete(limits = rev(levels(W$Feature)))

# as.character(W$Feature)


