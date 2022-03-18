#!/usr/bin/env Rscript

library(stringr)
library(parallel)

con = file("stdin", "r")
while(length(lines <- readLines(con, n = 100000)) > 0) { # 100k is best
    cat(paste0(mclapply(lines, function(line) {
        a = sort(table(str_split(line, "\t")), decreasing=TRUE)
        sprintf("ChromHMM;%s\t%s", names(a)[names(a)!="NA"][1],
            paste0(sprintf("%s;%d", names(a), a), collapse=","))
    }, mc.cores=24), collapse="\n"),"\n", sep="")
}
