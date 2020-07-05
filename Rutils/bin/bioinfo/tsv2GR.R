#!/usr/bin/env Rscript
library(stringr)
library(GenomicRanges)
args <- commandArgs(trailingOnly=TRUE)

tsvfile <- args[1]
grfile <- args[2]

## bed <- read.table(gzfile(tsvfile), sep='\t', header=T, row.names='probeID')
bed <- read.table(gzfile(tsvfile), sep='\t', header=T)
chrms <- c(paste0('chr',1:22),'chrX','chrY','chrM','*')
bed[!is.na(bed$beg),'beg'] <- bed[!is.na(bed$beg),'beg']+1
bed[is.na(bed$beg),'beg'] <- 0
bed[is.na(bed$end),'end'] <- 0
bed[is.na(bed$chrm),'chrm'] <- '*'
gr <- GRanges(seqnames=bed$chrm, ranges=IRanges(bed$beg, bed$end), strand=bed$probe_strand, seqinfo=Seqinfo(chrms))
mcols(gr) <- bed[,!(colnames(bed) %in% c('chrm','beg','end','probe_strand','cpg'))]
names(gr) <- bed$cg
gr <- sort(gr, by = ~ seqnames + start + end)
saveRDS(gr, file=grfile) # compress='xz' is about half the size compared to gz, but they are not directly downloadable through gzcon(url())
