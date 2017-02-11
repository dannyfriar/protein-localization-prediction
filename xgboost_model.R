library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(xgboost)
library(caret)
library(scales)

#----------- Read data, split out validation set and split out labels
train <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')
test <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')

# # Use important features only
# train <- train[, c('isolectric_point', 'E_count', 'K_count', 'E_first50_count',  'aromaticity','L_first50_count',
#                    'molecular_weight', 'C_count','L_count', 'G_count', 'P_count', 'Q_count', 'D_first50_count',
#                    'N_count', 'I_count', 'M_count','K_first50_count',  'sequence_length',  'G_first50_count',
#                    'H_count', 'A_count', 'D_count', 'A_first50_count', 'C_first50_count', 'F_count', 'name', 'sequence', 'class')]
# 
# test <- test[, c('isolectric_point', 'E_count', 'K_count', 'E_first50_count',  'aromaticity','L_first50_count',
#                    'molecular_weight', 'C_count','L_count', 'G_count', 'P_count', 'Q_count', 'D_first50_count',
#                    'N_count', 'I_count', 'M_count','K_first50_count',  'sequence_length',  'G_first50_count',
#                    'H_count', 'A_count', 'D_count', 'A_first50_count', 'C_first50_count', 'F_count', 'name', 'sequence', 'class')]

smp_size <- floor(0.75 * nrow(train))
train_ind <- sample(seq_len(nrow(train)), size = smp_size)
train_set <- train[train_ind, ]
val_set <- train[-train_ind, ]

y_train <- as.numeric(factor(train_set$class, level=c('cyto', 'mito', 'nuclear', 'secreted'), labels=c(0, 1, 2, 3))) - 1
names_train <- train_set$name
X_train <- subset(train_set, select=-c(name, class, sequence))

y_val <- as.numeric(factor(val_set$class, level=c('cyto', 'mito', 'nuclear', 'secreted'), labels=c(0, 1, 2, 3))) - 1
names_val <- val_set$name
X_val <- subset(val_set, select=-c(name, class, sequence))

y_test <- as.numeric(factor(test$class, level=c('cyto', 'mito', 'nuclear', 'secreted'), labels=c(0, 1, 2, 3))) - 1
names_test <- test$name
X_test <- subset(test, select=-c(name, class, sequence))


#----------- XGBoost model
dtrain = xgb.DMatrix(as.matrix(X_train), label= y_train)
dval = xgb.DMatrix(as.matrix(X_val), label= y_val)

# XGBoost parameters
xgb_params = list(
  seed = 0,
  colsample_bytree = 0.7,
  subsample = 0.7,
  objective = 'multi:softmax',
  max_depth = 10,
  alpha = 1,
  lambda = 5,
  gamma = 0,
  min_child_weight = 20,
  num_class = 4
)

best_n_rounds = 100

# train data
gb_dt = xgb.train(xgb_params, dtrain, nrounds = as.integer(best_n_rounds))
train_pred = predict(gb_dt, dtrain)
print("Training Accuracy:")
print(sum(train_pred==y_train)/length(y_train))

# Test on validation set
val_pred = predict(gb_dt, dval)
print("Validation Accuracy:")
print(sum(val_pred==y_val)/length(y_val))

# Feature importance
importance_model <- xgb.importance(feature_names = names(X_train), model = gb_dt)
important_features <- importance_model$Feature[1:25]

importance_model$Feature <- factor(importance_model$Feature)
importance_model$Feature <- reorder(importance_model$Feature, -importance_model$Gain)

g <- ggplot(data=importance_model, aes(x=Feature, y=Gain, fill=Feature)) + geom_bar(stat='identity')
g <- g + labs(x='Feature', y='Contribution to Model') + coord_flip() + guides(fill=FALSE)
g + scale_x_discrete(limits = rev(levels(importance_model$Feature))) + scale_y_continuous(breaks=seq(0, 0.1, 0.01), labels=percent)