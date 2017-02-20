library(rvest)
library(plyr)
library(dplyr)
library(data.table)

url <- read_html('http://genome.unmc.edu/LocSigDB/cgi-bin/search.cgi?signal=all')
df <- url %>% html_node('table') %>% html_table(fill=TRUE)

nls_signals <- subset(data.table(df), select=c(X1, X3))
names(nls_signals) <- c('signal', 'location')
nls_signals <- nls_signals[!grep(' ', nls_signals$signal)]
nls_signals <- nls_signals[!grep('}', nls_signals$signal)]

nls_signals_vector <- nls_signals$signal
write.csv(nls_signals_vector, file = "~/Desktop/CSML/bioinformatics/coursework/data/NLS_signals.csv", row.names=FALSE)