#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
system("zcat %s | awk '!/^!/ && length($0)>0' >1", args[1])
## zcat GSE43975_series_matrix.txt.gz | awk '/^!Sample/' | wzmanip transpose - | cut -f1,2,8,10-13 >samplesheet.tsv

