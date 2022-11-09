#!/usr/bin/env Rscript
# suppressMessages(library(tidyverse))
suppressMessages(library(readr))
suppressMessages(library(readxl))
args <- commandArgs(trailingOnly=TRUE)

# '/Users/zhouw3/repo/samplesheets/201x/2016_04_05_TCGA_merged_mapping.xlsx'
if (length(args)>1){
	df <- read_excel(args[1], args[2])
} else {
	df <- read_excel(args[1])
}
cat(format_tsv(df))
