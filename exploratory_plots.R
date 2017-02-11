library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(xgboost)
library(caret)
library(scales)

#----------- Read data, split out validation set and split out labels
train <- read.csv('~/Desktop/CSML/bioinformatics/coursework/data/train.csv')
train$class <- factor(train$class, levels = c('cyto', 'secreted', 'mito', 'nuclear'), labels = c('Cytoplasm', 'Secreted', 'Mitochondria', 'Nuclear'))

# Plot of class
g <- ggplot(data=train, aes(x=class, fill=class)) + geom_bar() + labs(x='Subcelluar Location', y='Count')
g + guides(fill=FALSE) + scale_y_continuous(breaks=seq(0, 3000, 500)) 

# Plot of class against molecular weights
g <- ggplot(data=train, aes(x=class, y=molecular_weight, fill=class)) + geom_boxplot() + labs(y='Molecular Weight (Da)', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 500000)

# Plot of isoelectric point against molecular weights
g <- ggplot(data=train, aes(x=class, y=isolectric_point, fill=class)) + geom_boxplot() + labs(y='Isolectric Point', x='Subcelluar Location') + guides(fill=FALSE)
g

# Plot of aromaticity against molecular weights
g <- ggplot(data=train, aes(x=class, y=aromaticity, fill=class)) + geom_boxplot() + labs(y='Aromaticity', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 0.3)

# Plot of E count against molecular weights
g <- ggplot(data=train, aes(x=class, y=E_count, fill=class)) + geom_boxplot() + labs(y='E Count', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 150)

# Plot of K count against molecular weights
g <- ggplot(data=train, aes(x=class, y=K_count, fill=class)) + geom_boxplot() + labs(y='K Count', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 150)

# Heat map of last amino acid
train$sequence <- as.character(train$sequence)
train$last_aa <- substr(train$sequence, nchar(train$sequence), nchar(train$sequence))
last_class <- data.table(subset(train, select=c(class, last_aa)))
last_class$last_aa <- factor(last_class$last_aa)
last_class$count <- 1
last_class <- last_class[ , lapply(.SD, sum),  by = c('class', 'last_aa')][order(class, -count), ]
last_class[class=='Cytoplasm']$count <- last_class[class=='Cytoplasm']$count / sum(last_class[class=='Cytoplasm']$count)
last_class[class=='Secreted']$count <- last_class[class=='Secreted']$count / sum(last_class[class=='Secreted']$count)
last_class[class=='Mitochondria']$count <- last_class[class=='Mitochondria']$count / sum(last_class[class=='Mitochondria']$count)
last_class[class=='Nuclear']$count <- last_class[class=='Nuclear']$count / sum(last_class[class=='Nuclear']$count)
ggplot(last_class, aes(class, last_aa)) + geom_tile(aes(fill = count), colour = "white") + scale_fill_gradient(low = "white", high = "red")

