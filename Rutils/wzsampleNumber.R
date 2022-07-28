#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
n = as.integer(args[1])
m = as.integer(args[2])

sampled = sort(sample(n, m, replace=FALSE))
cat(paste0(sampled, collapse='\n'),'\n')
