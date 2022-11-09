#!/usr/bin/env Rscript

library(stringr)
library(parallel)

con = file("stdin", "r")
while(length(lines <- readLines(con, n = 100000)) > 0) { # 100k is best
a = do.call(rbind, lapply(lines, function(x) as.numeric(str_split(x, "\t")[[1]])))
cat(paste0(paste0(rowSums(a == 1, na.rm=T), "\t", rowSums(a == 0)), collapse="\n"))
}

