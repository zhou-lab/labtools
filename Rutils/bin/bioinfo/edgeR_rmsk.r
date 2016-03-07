#!/usr/bin/env r
## this is to compare counts from rmsk mapping from base-wise counting
## edge_rmsk.r -a condition1 -b condition2 -A condition1bam.rep1,condition2bam.rep2 -B condition2bam.rep1,condition2bam.rep2 -o edgeR_rmsk.output
suppressPackageStartupMessages(library(docopt))

"Usage:
  edgeR_rmsk.r -a CONDITION1 -b CONDITION2 -A BAMS1 -B BAMS2 -o OUTPUT

Options:
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -A BAMS1        path to bams for condition1 (comma-separated)
  -B BAMS2        path to bams for condition2 (comma-separated)
  -o OUTPUT       output file path
" -> doc

opt <- docopt(doc)
## cat(str(opt))

cat(format(Sys.time()), "Reading inputs ...\n")
bam.paths1 <- unlist(strsplit(opt$A,','))
bam.paths2 <- unlist(strsplit(opt$B,','))
n1 <- length(bam.paths1)
n2 <- length(bam.paths2)
t <- sapply(c(bam.paths1, bam.paths2), read.table, simplify=F)

suppressPackageStartupMessages(library(edgeR))

cat(format(Sys.time()), "Normalizing ...\n")
toc <- sapply(t, function(x) {ifelse(x$V4=="-", -x$V9, x$V8)})
normFactors = calcNormFactors(as.matrix(toc))
metacol <- t[[1]][,1:7]
row.names(toc) <- paste0(metacol$V1,"_",metacol$V2,"_",metacol$V3,"_",metacol$V5)

group <- rep(c(opt$a, opt$b), c(n1, n2))
DGE <- DGEList(toc, lib.size=normFactors*colSums(toc), group=group)

cat(format(Sys.time()), "Estimating dispersion ...\n")
if (n1>1 && n2>1) { # with replica
    DGE.withdisp=estimateCommonDisp(DGE)
} else {                                # no replica
    cat("\nWarning: no replica, use GLM to estimate dispersion...")
    DGE.withdisp=estimateGLMCommonDisp(DGE,robust=T,subset=NULL,method="deviance")
}
tested=exactTest(DGE.withdisp)

cat(format(Sys.time()), "Outputing ...\n")
write.table(cbind(metacol, DGE.withdisp$counts, tested$table), file=opt$o,sep="\t", quote=F, col.names=NA)

cat(format(Sys.time()), "Done.")

