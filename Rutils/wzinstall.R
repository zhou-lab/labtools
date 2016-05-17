#!/usr/bin/env Rscript
library(R6)
target <- commandArgs(trailingOnly=TRUE)

library(devtools)
install_local(target)
