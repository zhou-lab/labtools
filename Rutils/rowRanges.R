args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(dplyr))
suppressMessages(library(SummarizedExperiment))
options(scipen=99)

write.table(as.data.frame(rowRanges(readRDS(args[[1]]))),file="",sep="\t",quote=FALSE, row.names=FALSE)

