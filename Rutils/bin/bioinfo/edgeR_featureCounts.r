#!/usr/bin/env r
## edgeR_featureCounts.r -g mm10 -a condition1 -b condition2 -A condition1bam.rep1,condition2bam.rep2 -B condition2bam.rep1,condition2bam.rep2 -o edgeRoutput
suppressPackageStartupMessages(library(docopt))

"Usage:
  edgeR.r [-G GTF] [-P] [-F CNTS ] [-s STARTCOL] [-T NTHREADS] -a CONDITION1 -b CONDITION2 -A BAMS1 -B BAMS2 -o OUTPUT

Options:
  -G GTF          Ensembl GTF file (if not provided, downloaded from UCSC)
  -F CNTS         count table from FeatureCounts, column names must match those in -A and -B
  -s STARTCOL     starting column index after annotation [default: 7]
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

cat(format(Sys.time()), "Making table of counts...\n")
suppressPackageStartupMessages(library(Rsubread))
samples1 <- unlist(strsplit(opt$A,','))
samples2 <- unlist(strsplit(opt$B,','))
n1 <- length(samples1)
n2 <- length(samples2)

##### below performs feature counts in R
## NOTE: edgeR needs raw counts, not FPKM
if (!is.null(opt$G)) {
    cat(format(Sys.time()), "Run featureCounts...\n")
    cnts <- featureCounts(files=c(samples1, samples2), annot.inbuilt=opt$G, nthreads=as.numeric(opt$T), isPairedEnd=opt$P, annot.ext=opt$G, ignoreDup=T, isGTFAnnotationFile=T)
    toc <- cnts$counts
} else if (!is.null(opt$F)) {
    cnts <- read.table(opt$F,stringsAsFactors=F, header=T, check.names=F)
    toc <- cnts[,c(samples1,samples2)]
} else {
    stop("Either provide GTF file to work on bams (-G) or provide featureCounts file (-F)")
}
## write.table(merge(cnts$counts, cnts$annotation, by.x=0, by.y="GeneID"), file=opt$o, quote=F, sep="\t", row.names=F)

cat(format(Sys.time()), "Normalizing...\n")
suppressPackageStartupMessages(library(edgeR))
normFactors = calcNormFactors(as.matrix(toc))

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

cat(format(Sys.time()), "Outputing...\n")

if (is.null(opt$F)) {
    write.rowname <- TRUE
    write.colname <- NA
    meta <- DGE.withdisp$counts
} else {
    write.rowname <- FALSE
    write.colname <- TRUE
    meta <- cbind(cnts[,1:(as.numeric(opt$s)-1)], cnts[,c(samples1,samples2)])
}

out <- cbind(meta, tested$table)
out.ordered <- out[order(out$PValue),]
write.table(out.ordered, file=opt$o, sep="\t", quote=F, row.names=write.rowname, col.names=write.colname)

cat(format(Sys.time()), "Done.\n")
sink()
