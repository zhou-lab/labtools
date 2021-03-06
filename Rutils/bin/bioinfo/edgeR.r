#!/usr/bin/env r
## this is memory heavy, since it loads the bam into memory, allow 10G
## this is using GenomicAlignment, GenomicFeatures and GRangesList packages which is quite unwieldy
## no filtering of reads based on mapping quality
## EdgeR.r -g mm10 -a condition1 -b condition2 -A condition1bam.rep1,condition2bam.rep2 -B condition2bam.rep1,condition2bam.rep2 -o edgeRoutput
suppressMessages(library(docopt))

"Usage:
  edgeR.r -g REFVERSION [-G GTF] -a CONDITION1 -b CONDITION2 -A BAMS1 -B BAMS2 -o OUTPUT

Options:
  -g REFVERSION   reference version (e.g, hg19, mm10)
  -G GTF          Ensembl GTF file (if not provided, downloaded from UCSC)   
  -a CONDITION1   name of condition1
  -b CONDITION2   name of condition2
  -A BAMS1        path to bams for condition1 (comma-separated)
  -B BAMS2        path to bams for condition2 (comma-separated)
  -o OUTPUT       output file path
" -> doc

opt <- docopt(doc)
sink(stderr())
## cat(str(opt))

cat(format(Sys.time()), "Retrieving annotations...\n")
suppressPackageStartupMessages(library(GenomicFeatures))
if (is.null(opt$G)) {
    txdb <- suppressWarnings(makeTxDbFromUCSC(genome=opt$g,tablename='ensGene'))
    gtf <- NULL
} else {
    txdb <- makeTxDbFromGFF(opt$G, format="gtf",circ_seqs=character())
    seqlevelsStyle(txdb) <- "UCSC"
    suppressPackageStartupMessages(library(rtracklayer))
    gtf <- import(opt$G)
}

## txdb is TxDb class from GenomicFeatures package
## txByGene is GRangesList
txByGene <- transcriptsBy(txdb,'gene')

## NOTE: edgeR needs raw counts, not FPKM
cat(format(Sys.time()), "Making table of counts...\n")
suppressPackageStartupMessages(library(GenomicAlignments))
toc1 <- data.frame(sapply(unlist(strsplit(opt$A,',')), function(bam.path) {
    cat("Loading", bam.path, "..."); flush.console()
    gal <- readGAlignments(bam.path)
    reads <- GRanges(seqnames = rname(gal), ranges=IRanges(start=start(gal),end=end(gal)), strand=strand(gal))
    ## unstranded library
    ## reads=GRanges(seqnames=rname(reads),ranges=IRanges(start=start(reads),end=end(reads)), strand=rep("*",length(reads)))
    counts <- countOverlaps(txByGene, reads)
    cat("OK\n"); flush.console()
    counts
}, USE.NAMES=T), stringsAsFactors=F)

toc2 <- data.frame(sapply(unlist(strsplit(opt$B,',')), function(bam.path) {
    cat("Loading", bam.path, "..."); flush.console()
    gal <- readGAlignments(bam.path)    # GAlignments objects
    reads <- GRanges(seqnames = rname(gal), ranges=IRanges(start=start(gal),end=end(gal)), strand=strand(gal))
    ## unstranded library
    ## reads=GRanges(seqnames=rname(reads),ranges=IRanges(start=start(reads),end=end(reads)), strand=rep("*",length(reads)))
    counts <- countOverlaps(txByGene, reads)
    cat("OK\n"); flush.console()
    counts
}, USE.NAMES=T), stringsAsFactors=F)

toc <- cbind(toc1, toc2)

cat(format(Sys.time()), "Normalizing ...")
suppressPackageStartupMessages(library(edgeR))
normFactors = calcNormFactors(as.matrix(toc))
cat("OK\n")

cat(format(Sys.time()), "Statistical testing...\n")
group <- rep(c(opt$a, opt$b), c(dim(toc1)[2],dim(toc2)[2]))
DGE <- DGEList(toc, lib.size=normFactors*colSums(toc), group=group)
if (dim(toc1)[2]>1 && dim(toc2)[2]>1) { # with replica
    DGE.withdisp=estimateCommonDisp(DGE)
} else {                                # no replica
    cat("\nWarning: no replica, use GLM to estimate dispersion...")
    DGE.withdisp=estimateGLMCommonDisp(DGE,robust=T,subset=NULL,method="deviance")
}
tested=exactTest(DGE.withdisp)

cat(format(Sys.time()), "Outputing ...\n")
if (is.null(gtf)) {
    ## Old way of retrieving annotation from OrgDb, this lacks non-coding element names
    if (opt$g == 'mm10') {
        suppressPackageStartupMessages(library(org.Mm.eg.db))
        annodb <- org.Mm.eg.db
    } else if (opt$g == 'hg19') {
        suppressPackageStartupMessages(library(org.Hs.eg.db))
        annodb <- org.Hs.eg.db
    }
    symbol <- mapIds(annodb, keys=row.names(tested$table), column="SYMBOL",keytype="ENSEMBL",multiVals="first")
} else {
    meta <- unique(as.data.frame(elementMetadata(gtf)[,c("gene_id","gene_name")]))
    id2name <- new.env(parent=emptyenv())
    invisible(apply(meta,1,function(x) id2name[[x[1]]] <<- x[[2]]))
    symbol <- unlist(mget(row.names(tested$table), envir=id2name));
}
write.table(cbind(DGE.withdisp$counts, tested$table, data.frame(symbol=symbol)), file=opt$o,sep="\t", quote=F, col.names=NA)

cat(format(Sys.time()), "Done.")
sink()
