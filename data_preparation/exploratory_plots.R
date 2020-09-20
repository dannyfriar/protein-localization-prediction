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
# 
# train <- data.table(train)
# train <- subset(train[class=='Nuclear'], select=c(class, sequence))
# train$sequence <- as.character(train$sequence)
# train$end_seq <- substr(train$sequence, nchar(train$sequence)-9, nchar(train$sequence))
# train$first_seq <- substr(train$sequence, 1, 10)
# train$sequence <- NULL

# Plot of class
g <- ggplot(data=train, aes(x=class, fill=class)) + geom_bar() + labs(x='Subcelluar Location', y='Count')
g + guides(fill=FALSE) + scale_y_continuous(breaks=seq(0, 3000, 500))

# Plot of class against molecular weights
g <- ggplot(data=train, aes(x=class, y=molecular_weight, fill=class)) + geom_boxplot() + labs(y='Molecular Weight (Da)', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 200000)

# Plot of isoelectric point against molecular weights
g <- ggplot(data=train, aes(x=class, y=isolectric_point, fill=class)) + geom_boxplot() + labs(y='Isolectric Point', x='Subcelluar Location') + guides(fill=FALSE)
g

# Plot of pct_hydrophobic against molecular weights
g <- ggplot(data=train, aes(x=class, y=pct_hydrophobic, fill=class)) + geom_boxplot() + labs(y='Percent Hydrophobic', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0.45, 0.75)

# Plot NLS signals counts against class
g <- ggplot(data=train, aes(x=class, y=nls_count, fill=class)) + geom_boxplot() + labs(y='NLS signals', x='Subcelluar Location') + guides(fill=FALSE)
g

# Plot of sequence length against class
g <- ggplot(data=train, aes(x=class, y=sequence_length, fill=class)) + geom_boxplot() + labs(y='Sequence Length', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 1500)

# Plot of aromaticity against molecular weights
g <- ggplot(data=train, aes(x=class, y=aromaticity, fill=class)) + geom_boxplot() + labs(y='Aromaticity', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 0.3)

# Plot of C count against molecular weights
g <- ggplot(data=train, aes(x=class, y=C_count, fill=class)) + geom_boxplot() + labs(y='Proportion of Cysteine AAs', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 0.06)

# Plot of K count against molecular weights
g <- ggplot(data=train, aes(x=class, y=K_count, fill=class)) + geom_boxplot() + labs(y='K Count', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 0.15)

# Percent +ve charge
g <- ggplot(data=train, aes(x=class, y=pct_pos_charged, fill=class)) + geom_boxplot() + labs(y='Pct +ve charged', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 0.25)

# Percent -ve charge
g <- ggplot(data=train, aes(x=class, y=pct_neg_charged, fill=class)) + geom_boxplot() + labs(y='Pct -ve charged', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0, 0.25)

# Percent hydrophobic
g <- ggplot(data=train, aes(x=class, y=pct_hydrophobic, fill=class)) + geom_boxplot() + labs(y='Pct hydrophobic', x='Subcelluar Location') + guides(fill=FALSE)
g + ylim(0.25, 0.8)

# # Heat map of last amino acids
# train$sequence <- as.character(train$sequence)
# train$last_aa <- substr(train$sequence, nchar(train$sequence)-1, nchar(train$sequence))
# last_class <- data.table(subset(train, select=c(class, last_aa)))
# last_class$last_aa <- factor(last_class$last_aa)
# last_class$count <- 1
# last_class <- last_class[ , lapply(.SD, sum),  by = c('class', 'last_aa')][order(class, -count), ]
# last_class[class=='Cytoplasm']$count <- last_class[class=='Cytoplasm']$count / sum(last_class[class=='Cytoplasm']$count)
# last_class[class=='Secreted']$count <- last_class[class=='Secreted']$count / sum(last_class[class=='Secreted']$count)
# last_class[class=='Mitochondria']$count <- last_class[class=='Mitochondria']$count / sum(last_class[class=='Mitochondria']$count)
# last_class[class=='Nuclear']$count <- last_class[class=='Nuclear']$count / sum(last_class[class=='Nuclear']$count)
# last_class <- last_class[count >=0.01]
# ggplot(last_class, aes(class, last_aa)) + geom_tile(aes(fill = count), colour = "white") + scale_fill_gradient(low = "white", high = "red")
# 
# 
# # Heat map of first amino acids
# train$sequence <- as.character(train$sequence)
# train$first_aa <- substr(train$sequence, 2, 2)
# first_class <- data.table(subset(train, select=c(class, first_aa)))
# first_class$first_aa <- factor(first_class$first_aa)
# first_class$count <- 1
# first_class <- first_class[ , lapply(.SD, sum),  by = c('class', 'first_aa')][order(class, -count), ]
# first_class[class=='Cytoplasm']$count <- first_class[class=='Cytoplasm']$count / sum(first_class[class=='Cytoplasm']$count)
# first_class[class=='Secreted']$count <- first_class[class=='Secreted']$count / sum(first_class[class=='Secreted']$count)
# first_class[class=='Mitochondria']$count <- first_class[class=='Mitochondria']$count / sum(first_class[class=='Mitochondria']$count)
# first_class[class=='Nuclear']$count <- first_class[class=='Nuclear']$count / sum(first_class[class=='Nuclear']$count)
# first_class <- first_class[count >=0.01]
# ggplot(first_class, aes(class, first_aa)) + geom_tile(aes(fill = count), colour = "white") + scale_fill_gradient(low = "white", high = "red")

