#!/usr/bin/env Rscript-4.3.1
## 202307_SE2bed.R "DNAm/SE/{arg}.rds"
## ld22c()
suppressMessages(library(readr))
suppressMessages(library(SummarizedExperiment))

args = commandArgs(trailingOnly=TRUE)
se = readRDS(args[1])
df = as.data.frame(rowData(se))
df = cbind(chrm=df$chrm, beg=df$beg, end=df$beg+2, df[,3:ncol(df)])
cat(format_tsv(df))

## cd = as.data.frame(read_tsv("DNAm/pseudos.cg.idx", col_names=c("celltype", "addr")))
## mtx = as.matrix(read.table(text=system(sprintf("kycg unpack -a %s", args[1]), intern=TRUE), header=F))
## colnames(mtx) = rownames(cd)
## gr = bed2GR(args[2], rangeOnly=T)
## names(gr) = paste0(seqnames(gr),"_",start(gr))
## se = SummarizedExperiment(mtx, colData=cd, rowRanges=gr)
## saveRDS(se, file=args[3])

