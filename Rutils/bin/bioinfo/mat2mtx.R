#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 2)
library(rmatio)
library(Matrix) ## sparse matrix dgCMatrix
m <- read.mat(args[[1]])
writeMM(m, file=paste(args[[2]]))
