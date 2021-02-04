#!/usr/bin/env Rscript
library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)

df <- readRDS(args[1])
write_tsv(df, args[2])
