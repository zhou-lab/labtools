#!/usr/bin/env Rscript-4.2.devel

target <- commandArgs(trailingOnly=TRUE)
format <- "html"
if (length(target) > 1) {
  format <- target[2]
  target <- target[1]
}

if (dir.exists(target)) {
  library(devtools)
  document()
}
