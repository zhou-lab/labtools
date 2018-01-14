#!/usr/bin/env Rscript
library(stringr)
library(GenomicRanges)
args <- commandArgs(trailingOnly=TRUE)

tsvfile <- args[1]
grfile <- args[2]

bed <- read.table(gzfile(tsvfile), sep='\t', header=T, row.names='probeID')
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM','*')
bed[!is.na(bed$CpG_beg),'CpG_beg'] <- bed[!is.na(bed$CpG_beg),'CpG_beg']+1
bed[is.na(bed$CpG_beg),'CpG_beg'] <- 0
bed[is.na(bed$CpG_end),'CpG_end'] <- 0
bed[is.na(bed$CpG_chrm),'CpG_chrm'] <- '*'
gr <- GRanges(seqnames=bed$CpG_chrm, ranges=IRanges(bed$CpG_beg, bed$CpG_end), strand=bed$probe_strand, seqinfo=Seqinfo(chrms))
mcols(gr) <- bed[,!(colnames(bed) %in% c('CpG_chrm','CpG_beg','CpG_end','probe_strand'))]
names(gr) <- rownames(bed)
gr <- sort(gr, by = ~ seqnames + start + end)
saveRDS(gr, file=grfile, compress='xz')
