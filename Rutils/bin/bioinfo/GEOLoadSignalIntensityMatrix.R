#!/usr/bin/env Rscript
library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)

## df <- read_csv(gzfile('GSE43414_signal_intensities.csv.gz'))
## a <- table(sapply(strsplit(colnames(df),' '), function(x) x[1]))
## samples <- names(a[a == 3])
## betas <- do.call(cbind, lapply(samples[1:3], function(x) {
##     df[[paste0(x, ' Methylated Signal')]] / (
##         df[[paste0(x, ' Unmethylated Signal')]] + df[[paste0(x, ' Methylated Signal')]])
## }))

## GEOLoad_SignalIntensityMatrix.R /mnt/isilon/zhoulab/20191212_GEO_datasets/GSE43414_multitissue/GSE43414_daten1_geo_all_cohorts.csv.gz
df <- read_csv(gzfile(args[[1]]))
betas <- do.call(cbind, lapply(seq(2, ncol(df), by = 2), function(i) df[[i]]))
pvals <- do.call(cbind, lapply(seq(3, ncol(df), by = 2), function(i) df[[i]]))
samples <- colnames(df)[seq(2, ncol(df), by = 2)]
rownames(betas) <- df[[1]]
rownames(pvals) <- df[[1]]
colnames(betas) <- samples
colnames(pvals) <- samples
saveRDS(betas, file='output_betas.rds')
saveRDS(pvals, file='output_pvals.rds')
