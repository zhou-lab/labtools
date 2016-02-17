#!/usr/bin/env r
## EdgeR_featureCounts.r -g mm10 -a condition1 -b condition2 -A condition1bam.rep1,condition2bam.rep2 -B condition2bam.rep1,condition2bam.rep2 -o edgeRoutput
suppressMessages(library(docopt))

"Usage:
  edgeR.r -G GTF [-P] [-T NTHREADS] -a CONDITION1 -b CONDITION2 -A BAMS1 -B BAMS2 -o OUTPUT

Options:
  -G GTF          Ensembl GTF file (if not provided, downloaded from UCSC)   
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -T NTHREADS     number of threads [default: 3]
  -A BAMS1        path to bams for condition1 (comma-separated)
  -B BAMS2        path to bams for condition2 (comma-separated)
  -o OUTPUT       output file path
  -P              whether bam is paired-end
" -> doc

opt <- docopt(doc)
sink(stderr())
## cat(str(opt))

## NOTE: edgeR needs raw counts, not FPKM
suppressPackageStartupMessages(library(Rsubread))
cat(format(Sys.time()), "Making table of counts...\n")
samples1 <- unlist(strsplit(opt$A,','))
samples2 <- unlist(strsplit(opt$B,','))
n1 <- length(samples1)
n2 <- length(samples2)
cnts <- featureCounts(files=c(samples1, samples2), annot.inbuilt=opt$G, nthreads=as.numeric(opt$T), isPairedEnd=opt$P, annot.ext=opt$G, ignoreDup=T, isGTFAnnotationFile=T)
toc <- cnts$counts

# write.table(merge(cnts$counts, cnts$annotation, by.x=0, by.y="GeneID"), file=opt$o, quote=F, sep="\t", row.names=F)

cat(format(Sys.time()), "Normalizing ...")
suppressPackageStartupMessages(library(edgeR))
normFactors = calcNormFactors(as.matrix(toc))
cat("OK\n")

cat(format(Sys.time()), "Statistical testing...\n")
group <- rep(c(opt$a, opt$b), c(n1, n2))
DGE <- DGEList(toc, lib.size=normFactors*colSums(toc), group=group)
if (n1>1 && n2>1) { # with replica
    DGE.withdisp=estimateCommonDisp(DGE)
} else {                                # no replica
    cat("\nWarning: no replica, use GLM to estimate dispersion...")
    DGE.withdisp=estimateGLMCommonDisp(DGE,robust=T,subset=NULL,method="deviance")
}
tested=exactTest(DGE.withdisp)

cat(format(Sys.time()), "Outputing ...\n")
write.table(cbind(DGE.withdisp$counts, tested$table), file=opt$o,sep="\t", quote=F, col.names=NA)

cat(format(Sys.time()), "Done.")
sink()
