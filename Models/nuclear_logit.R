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

# Load data
train <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')
load(file = "~/Desktop/CSML/bioinformatics/coursework/important_features.RData")
train <- train[, important_features]
train$class_binary <- 0
train <- data.table(train)
train[class=='nuclear']$class_binary <- 1
train <- subset(train, select=-c(class, name, sequence))
train$class_binary <- factor(train$class_binary)

# Logistic regression
glm(formula = class_binary ~ ., family = binomial(link = "logit"), data = train)
