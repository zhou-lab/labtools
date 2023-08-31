#!/usr/bin/env Rscript-4.3.1
## 202306_cgToSE <.cg file> <rowdata.bed> <output.rds>

ld22c()

args = commandArgs(trailingOnly=TRUE)
cd = as.data.frame(read_tsv(sprintf("%s.idx", args[1]), col_names=c("celltype", "addr")))
mtx = as.matrix(read.table(text=system(sprintf("kycg unpack -a %s", args[1]), intern=TRUE), header=F))
colnames(mtx) = rownames(cd)
gr = bed2GR(args[2], rangeOnly=T)
names(gr) = paste0(seqnames(gr),"_",start(gr))
se = SummarizedExperiment(mtx, colData=cd, rowRanges=gr)
saveRDS(se, file=args[3])
