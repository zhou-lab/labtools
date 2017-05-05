#!/usr/bin/env Rscript

a <- read.table('samples_curated.txt', header=T, sep='\t', stringsAsFactors=FALSE)
load('betas.rda')
rownames(a) <- rownames(samples)
samples <- a
save(betas, samples, file='betas.rda')
