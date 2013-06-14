#!/usr/bin/env Rscript

argv <- commandArgs(T)

Data <- read.table(argv[1], sep = '\t', header = T, as.is = T)
gtCols <- grep(".GT", colnames(Data))
normalGT <- Data[,gtCols[1]]
tumorGT <- Data[,gtCols[2]]

write.table(subset(Data, normalGT == "0/0" & tumorGT == "0/1"), quote = F, sep = '\t', row.names = F)

