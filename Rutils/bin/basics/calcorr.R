#!/usr/bin/env Rscript
library(data.table)
dt <- fread('file:///dev/stdin')

cat(cor(dt[,1], dt[,2], use='na.or.complete', method='spearman'),'\n')
