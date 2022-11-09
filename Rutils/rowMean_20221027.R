#!/usr/bin/env Rscript

library(stringr)
library(parallel)

con = file("stdin", "r")
while(length(lines <- readLines(con, n = 100000)) > 0) { # 100k is best
a = do.call(rbind, lapply(lines, function(x) as.numeric(str_split(x, "\t")[[1]])))
a[a==2] <- NA
cat(paste0(paste0(rowMeans(a, na.rm=T), "\t", rowSums(!is.na(a))), collapse="\n"))
}

